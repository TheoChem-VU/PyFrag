from __future__ import unicode_literals

import os
import copy

from os.path import join as opj

from ..core.basejob import SingleJob
from ..core.results import Results
from ..core.settings import Settings
from ..core.errors import PlamsError
from ..tools.kftools import KFFile

__all__ = ['ADFJob', 'ADFResults', 'BANDJob', 'BANDResults', 'DFTBJob', 'DFTBResults', 'FCFJob', 'FCFResults', 'DensfJob', 'DensfResults']


#===================================================================================================
#===================================================================================================


class SCMResults(Results):
    """Abstract class gathering common mechanisms for results of all ADF Suite binaries."""
    _kfext = ''


    def collect(self):
        """Collect files present in the job folder. Use parent method from |Results|, then create an instance of |KFFile| for the main KF file and store it as ``_kf`` attribute.
        """
        Results.collect(self)
        kfname = self.job.name + self.__class__._kfext
        if kfname in self.files:
            self._kf = KFFile(opj(self.job.path, kfname))


    def refresh(self):
        """Refresh the contents of ``files`` list. Use parent method from |Results|, then look at all attributes that are instances of |KFFile| and check if they point to existing files. If not, try to reinstantiate them with current job path (that can happen while loading a pickled job after the entire job folder was moved).
        """
        Results.refresh(self)
        to_remove = []
        for attr,val in self.__dict__.items():
            if isinstance(val, KFFile) and os.path.dirname(val.path) != self.job.path:
                guessnewpath = opj(self.job.path, os.path.basename(val.path))
                if os.path.isfile(guessnewpath):
                    self.__dict__[attr] = KFFile(guessnewpath)
                else:
                    to_remove.append(attr)
        for i in to_remove:
            del self.__dict__[i]


    def readkf(self, section, variable):
        """readkf(section, variable)
        Read data from *section*/*variable* of the main KF file.

        The type of the returned value depends on the type of *variable* defined inside KF file. It can be: single int, list of ints, single float, list of floats, single boolean, list of booleans or string. """
        if self._kf:
            return self._kf.read(section, variable)
        raise FileError('KFFile of %s not present in %s' % (self.job.name, self.job.path))


    def newkf(self, filename):
        """newkf(filename)
        Create new |KFFile| instance using file *filename* in the job folder.

        Example usage::

            >>> res = someadfjob.run()
            >>> tape13 = res.newkf('$JN.t13')
            >>> print(tape13.read('Geometry', 'xyz'))

        """
        self.refresh()
        filename = filename.replace('$JN', self.job.name)
        if filename in self.files:
            return KFFile(opj(self.job.path, filename))
        else:
            raise FileError('File %s not present in %s' % (filename, self.job.path))


    def get_molecule(self, section, variable, unit='bohr', internal=False, n=1):
        """get_molecule(section, variable, unit='bohr', internal=False, n=1)
        Read molecule coordinates from *section*/*variable* of the main KF file.

        Returned |Molecule| instance is created by copying a molecule from associated |SCMJob| instance and updating atomic coordinates with values read from *section*/*variable*. The format in which coordinates are stored is not consistent for all programs or even for different sections of the same KF file. Sometimes coordinates are stored in bohr, sometimes in angstrom. The order of atoms can be either input order or internal order. These settings can be adjusted with *unit* and *internal* parameters. Some variables store more than one geometry, in those cases *n* can be used to choose the preferred one.
        """
        if self.job.molecule:
            m = self.job.molecule.copy()
            natoms = len(m)
            coords = self.readkf(section, variable)
            coords = [coords[i:i+3] for i in range(0,len(coords),3)]
            if len(coords) > natoms:
                coords = coords[(n-1)*natoms:n*natoms]
            if internal:
                mapping = self._int2inp()
                coords = [coords[mapping[i]-1] for i in range(len(coords))]
            for at,coord in zip(m,coords):
                at.move_to(coord, unit)
            return m
        raise PlamsError('get_molecule() failed. Corresponding job has no molecule associated')


    def _kfpath(self):
        """_kfpath()
        Return an absolute path to the main KF file.
        """
        return opj(self.job.path, self.job.name + self.__class__._kfext)


    def _settings_reduce(self):
        """_settings_reduce()
        When this object is present as a value in some |Settings| instance and string representation is needed, use the absolute path to the main KF file. See :meth:`Settings.__reduce__<scm.plams.settings.Settings.__reduce__>` for details.
        """
        return self._kfpath()


    def _export_attribute(self, attr, other):
        """_export_attribute(attr, other)
        If *attr* is a KF file take care of a proper path. Otherwise use parent method. See :meth:`Results._copy_to<scm.plams.results.Results._copy_to>` for details.
        """
        if isinstance(attr, KFFile):
            oldname = os.path.basename(attr.path)
            newname = Results._replace_job_name(oldname, self.job.name, other.job.name)
            newpath = opj(other.job.path, newname)
            return KFFile(newpath) if os.path.isfile(newpath) else None
        else:
            return Results._export_attribute(self, attr, other)



    def _int2inp(self):
        """_int2inp()
        Obtain mapping from internal atom order to the input one. Abstract method.
        """
        raise PlamsError('Trying to run an abstract method SCMResults._int2inp()')



class SCMJob(SingleJob):
    """Abstract class gathering common mechanisms for jobs with all ADF Suite binaries."""
    _result_type = SCMResults
    _top = ['title','units','define']
    _command = ''
    _subblock_end = 'subend'



    def get_input(self):
        """Transform all contents of ``setting.input`` branch into string with blocks, keys and values.

        On the highest level alphabetic order of iteration is modified: keys occuring in attribute ``_top`` are printed first.

        Automatic handling of ``molecule`` can be disabled with ``settings.ignore_molecule = True``.
        """

        def parse(key, value, indent=''):
            ret = ''
            if isinstance(value, Settings):
                ret += indent + key
                if '_h' in value:
                    ret += ' ' + value['_h']
                ret += '\n'

                i = 1
                while ('_'+str(i)) in value:
                    ret += parse('', value['_'+str(i)], indent+'  ')
                    i += 1

                for el in value:
                    if not el.startswith('_'):
                        ret += parse(el, value[el], indent+'  ')

                if indent == '':
                    ret += 'end\n'
                else:
                    ret += indent + self._subblock_end + '\n'
            elif isinstance(value, list):
                for el in value:
                    ret += parse(key, el, indent)
            elif isinstance(value, (SCMJob, SCMResults, KFFile)):
                ret += parse(key, value._settings_reduce(), indent)
            elif value is '' or value is True:
                ret += indent + key + '\n'
            else:
                value = str(value)
                ret += indent + key
                if key != '' and not value.startswith('='):
                    ret += ' '
                ret += value + '\n'
            return ret

        inp = ''
        use_molecule = ('ignore_molecule' not in self.settings) or (self.settings.ignore_molecule == False)
        if use_molecule:
            self._parsemol()
        for item in self._top:
            item = self.settings.input.find_case(item)
            if item in self.settings.input:
                inp += parse(item, self.settings.input[item]) + '\n'
        for item in self.settings.input:
            if item.lower() not in self._top:
                inp += parse(item, self.settings.input[item]) + '\n'
        inp += 'end input\n'
        if use_molecule:
            self._removemol()
        return inp


    def get_runscript(self):
        """Generate a runscript. Returned string is of the form::

            $ADFBIN/name [-n nproc] <jobname.in [>jobname.out]

        ``name`` is taken from the class attribute ``_command``. ``-n`` flag is added if ``settings.runscript.nproc`` exists. ``[>jobname.out]`` is used based on ``settings.runscript.stdout_redirect``.
        """
        s = self.settings.runscript
        ret = '$ADFBIN/'+self._command
        if 'nproc' in s:
            ret += ' -n ' + str(s.nproc)
        ret += ' <'+self._filename('inp')
        if s.stdout_redirect:
            ret += ' >'+self._filename('out')
        ret += '\n\n'
        return ret


    def check(self):
        """Check if ``termination status`` variable from ``General`` section of main KF file equals ``NORMAL TERMINATION``."""
        try:
            ret = (self.results.readkf('General', 'termination status').strip(' ') == 'NORMAL TERMINATION')
        except:
            return False
        return ret


    def _parsemol(self):
        """Process |Molecule| instance stored in ``molecule`` attribute and add it as relevant entries of ``settings.input`` branch. Abstract method."""
        raise PlamsError('Trying to run an abstract method SCMJob._parsemol()')


    def _removemol(self):
        """Remove from ``settings.input`` all entries added by :meth:`_parsemol`. Abstract method."""
        raise PlamsError('Trying to run an abstract method SCMJob._removemol()')


    def _settings_reduce(self):
        """When this object is present as a value in some |Settings| instance and string representation is needed, use the absolute path to the main KF file. See :meth:`Settings.__reduce__<scm.plams.settings.Settings.__reduce__>` for details."""
        return self.results._kfpath()

    @staticmethod
    def _atom_symbol(atom):
        """Return the atomic symbol of *atom*. Ensure proper formatting for ADFSuite input taking into account ``ghost`` and ``name`` attributes of *atom*."""
        smb = atom.symbol if atom.atnum > 0 else ''  #Dummy atom should have '' instead of 'Xx'
        if hasattr(atom, 'ghost') and atom.ghost:
            smb = ('Gh.'+smb).rstrip('.')
        if hasattr(atom, 'name'):
            smb = (smb+'.'+str(atom.name)).lstrip('.')
        return smb


#===================================================================================================
#===================================================================================================


class ADFResults(SCMResults):
    _kfext = '.t21'
    _rename_map = {'TAPE21':'$JN'+_kfext, 'TAPE13':'$JN.t13', 'TAPE10':'$JN.t10'}

    def _int2inp(self):
        aoi = self.readkf('Geometry', 'atom order index')
        n = len(aoi)//2
        return aoi[:n]


class ADFJob(SCMJob):
    _result_type = ADFResults
    _command = 'adf'

    def _parsemol(self):
        for i,atom in enumerate(self.molecule):
            smb = self._atom_symbol(atom)
            suffix = ''
            if hasattr(atom,'fragment'):
                suffix += 'f={fragment} '
            if hasattr(atom,'block'):
                suffix += 'b={block}'

            self.settings.input.atoms['_'+str(i+1)] = ('%5i'%(i+1)) + atom.str(symbol=smb, suffix=suffix)

    def _removemol(self):
        if 'atoms' in self.settings.input:
            del self.settings.input.atoms


#===================================================================================================
#===================================================================================================


class BANDResults(SCMResults):
    _kfext = '.runkf'
    _rename_map = {'RUNKF':'$JN'+_kfext}

    def _int2inp(self):
        return self.readkf('geometry', 'Atom map new order')


class BANDJob(SCMJob):
    _result_type = BANDResults
    _command = 'band'

    def _parsemol(self):
        for i,atom in enumerate(self.molecule):
            self.settings.input.atoms['_'+str(i+1)] = atom.str(symbol=self._atom_symbol(atom))

        if self.molecule.lattice:
            for i,vec in enumerate(self.molecule.lattice):
                self.settings.input.lattice['_'+str(i+1)] = '%16.10f %16.10f %16.10f'%vec

    def _removemol(self):
        if 'atoms' in self.settings.input:
            del self.settings.input.atoms
        if 'lattice' in self.settings.input:
            del self.settings.input.lattice


#===================================================================================================
#===================================================================================================


class DFTBResults(SCMResults):
    _kfext = '.rkf'
    _rename_map = {'dftb.rkf':'$JN'+_kfext}

    def _int2inp(self):
        return list(range(1, 1+len(self.job.molecule)))


class DFTBJob(SCMJob):
    _result_type = DFTBResults
    _command = 'dftb'
    _top = ['units', 'task']
    _subblock_end = 'end'

    def _parsemol(self):
        s = self.settings.input
        for i,atom in enumerate(self.molecule):
            s[s.find_case('system')]['atoms']['_'+str(i+1)] = atom.str(symbol=self._atom_symbol(atom))
        if self.molecule.lattice:
            for i,vec in enumerate(self.molecule.lattice):
                s[s.find_case('system')]['lattice']['_'+str(i+1)] = '%16.10f %16.10f %16.10f'%vec

    def _removemol(self):
        s = self.settings.input
        system = s.find_case('system')
        if system in s:
            if 'atoms' in s[system]:
                del s[system]['atoms']
            if 'lattice' in s[system]:
                del s[system]['lattice']



#===================================================================================================
#===================================================================================================


class FCFResults(SCMResults):
    _kfext = '.t61'
    _rename_map = {'TAPE61':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('FCFResults do not support get_molecule() method. You can get molecules from job1 or job2')


class FCFJob(SCMJob):
    """A class representing calculation of Franck-Condon factors using ``fcf`` program.

    Two new attributes are introduced: ``inputjob1`` and ``inputjob2``. They are used to supply KF files from previous runs to ``fcf`` program. The value can either be a string with a path to KF file or an instance of any type of |SCMJob| or |SCMResults| (in this case the path to corresponding KF file will be extracted automatically). If the value of ``inputjob1`` or ``inputjob2`` is ``None``, no automatic handling occurs and user needs to manually supply paths to input jobs using proper keywords placed in ``myjob.settings.input`` (``STATES`` or ``STATE1`` and ``STATE2``).


    The resulting ``TAPE61`` file is renamed to ``jobname.t61``.
    """
    _result_type = FCFResults
    _command = 'fcf'
    _top = ['states', 'state1', 'state2']

    def __init__(self, inputjob1=None, inputjob2=None, **kwargs):
        SCMJob.__init__(self,**kwargs)
        self.inputjob1 = inputjob1
        self.inputjob2 = inputjob2

    def _parsemol(self):
        if isinstance(self.inputjob1, str):
            self.settings.input.state1 = self.inputjob1
        elif hasattr(self.inputjob1, '_settings_reduce'):
            self.settings.input.state1 = self.inputjob1._settings_reduce()

        if isinstance(self.inputjob2, str):
            self.settings.input.state2 = self.inputjob2
        elif hasattr(self.inputjob2, '_settings_reduce'):
            self.settings.input.state2 = self.inputjob2._settings_reduce()

    def _removemol(self):
        if 'state1' in self.settings.input:
            del self.settings.input.state1
        if 'state2' in self.settings.input:
            del self.settings.input.state2


#===================================================================================================
#===================================================================================================


class DensfResults(SCMResults):
    _kfext = '.t41'
    _rename_map = {'TAPE41':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('DensfResults do not support get_molecule() method. You can get molecule from inputjob')


class DensfJob(SCMJob):
    """A class representing calculation of molecular properties on a grid using ``densf`` program.

    A new attribute ``inputjob`` is introduced to supply KF file from previously run job. The value can either be a string with a path to KF file or an instance of any type of |SCMJob| or |SCMResults| (in this case the path to corresponding KF file will be extracted automatically). If the value of ``inputjob`` is ``None``, no automatic handling occurs and user needs to manually supply path to input job using ``INPUTFILE`` keyword placed in ``myjob.settings.input``.

    The resulting ``TAPE41`` file is renamed to ``jobname.t41``.
    """
    _result_type = DensfResults
    _command = 'densf'
    _top = ['inputfile', 'units']

    def __init__(self, inputjob=None, **kwargs):
        SCMJob.__init__(self, **kwargs)
        self.inputjob = inputjob

    def _parsemol(self):
        if isinstance(self.inputjob, str):
            self.settings.input.inputfile = self.inputjob
        elif hasattr(self.inputjob, '_settings_reduce'):
            self.settings.input.inputfile = self.inputjob._settings_reduce()

    def _removemol(self):
        if 'inputjob' in self.settings.input:
            del self.settings.input.inputjob

    def check(self):
        try:
            grep = self.results.grep_file('$JN.err', 'NORMAL TERMINATION')
        except:
            return False
        return len(grep) > 0
