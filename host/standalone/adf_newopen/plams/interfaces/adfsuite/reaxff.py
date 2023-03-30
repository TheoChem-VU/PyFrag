import numpy as np
import os
import shutil

from os.path import join as opj

from ...core.basejob import SingleJob
from ...core.errors import FileError
from ...core.functions import log
from ...core.settings import Settings
from ...tools.units import Units
from .scmjob import SCMResults
from .scmjob import SCMJob, SCMResults


__all__ = ['ReaxFFJob', 'ReaxFFResults', 'load_reaxff_control', 'reaxff_control_to_settings']


class ReaxFFResults(SCMResults):
    _rename_map = {'reaxout.kf':'$JN.kf'}
    _kfext = '.kf'

    def _int2inp(self):
        return list(range(1, 1+len(self.job.molecule)))


class ReaxFFJob(SingleJob):
    _filenames = {'inp':'control', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}
    _result_type = ReaxFFResults
    ffield_path = opj('$AMSHOME','atomicdata','ForceFields','ReaxFF')
    default_cell_size = 100.0

    check = SCMJob.check


    def get_input(self):
        """Produce the ``control`` file based on key-value pairs present in ``settings.input.control`` branch."""

        s = self.settings.input.control
        ret = ''
        order = s._order if '_order' in s else iter(s)
        for key in order:
            ret += '{:>7} {:6}\n'.format(s[key], key)
        return ret


    def get_runscript(self):
        """Generate a runscript.

        Returned string is just ``$AMSBIN/reaxff``, possibly prefixed with ``export NSCM=(number)`` if ``settings.runscript.nproc`` is present.
        """
        s = self.settings.runscript
        ret = ''
        if 'nproc' in s:
            ret += 'export NSCM={}\n'.format(s.nproc)
        ret += '$AMSBIN/reaxff'
        if s.stdout_redirect:
            ret += ' >'+self._filename('out')
        ret += '\n\n'
        return ret


    def hash_input(self):
        """Disable hashing for ReaxFF jobs.

        It is a common task in molecular dynamics to run several trajectories with the same initial conditions. In such a case |RPM| would prevent second and all consecutive executions. Hence we decided to disable |RPM| for ReaxFF.

        If you wish to bring it back, simply put ``ReaxFFJob.hash_input = SingleJob.hash_inputs`` somehwere at the beginning of your script.
        """
        return None


    def _get_ready(self):
        """Prepare contents of the job folder for execution.

        Use the parent method from |SingleJob| to produce the runscript and the input file (``control``). Then create ``ffield`` and ``geo`` files using, respectively, :meth:`_write_ffield` and :meth:`_write_geofile`.

        Then copy to the job folder all files listed in ``settings.input.external``. The value of this key should either be a list of strings with paths to files or a dictionary (also |Settings|) with paths to files as values and names under which these files should be copied to the job folder as keys.
        """
        SingleJob._get_ready(self)
        self._write_ffield(self.settings.input.ffield)
        self._write_geofile(molecule=self.molecule, filename=opj(self.path, 'geo'), settings=self.settings.input.geo, description=self.name, lattice=True)

        if 'external' in self.settings.input:
            ext = self.settings.input.external
            if isinstance(ext, list):
                for val in ext:
                    if os.path.isfile(val):
                        shutil.copy(val, self.path)
                    else:
                        raise FileError('File {} not present'.format(val))
            elif isinstance(ext, dict):
                for key, val in ext.items():
                    if os.path.isfile(val):
                        shutil.copy(val, opj(self.path, key))
                    else:
                        raise FileError('File {} not present'.format(val))



    def _write_ffield(self, ffield):
        """Copy to the job folder a force field file indicated by *ffield**.

        *ffield* should be a string with a path to some external file or with a filename present in ``$AMSHOME/atomicdata/ForceFields/ReaxFF``. The location of this search folder is defined by ``ffield_path`` class attribute).

        Given file is always coied to the job folder as ``ffield``, due to ReaxFF program requirements.
        """

        if os.path.isfile(ffield):
            shutil.copy(ffield, opj(self.path, 'ffield'))
        else:
            path = os.path.expandvars(opj(self.ffield_path, ffield))
            if os.path.isfile(path):
                shutil.copy(path, opj(self.path, 'ffield'))
            else:
                raise FileError('settings.input.ffield={} is neither a path to a file nor an existing force field from {}'.format(ffield, self.ffield_path))



    def _write_geofile(self, molecule, filename, settings, description, lattice=False):
        """Write to *filename* a geo-file describing *molecule*.

        *settings* should be a |Settings| instance containing all the additional key-value pairs that should be present in the resulting geo-file. To obtain multiple occurrences of the same key in the geo-file, put all the values as a list in *settings*.

        *description* is the default value for ``DESCRP`` key. It is used only if ``descrp`` key is not present in *settings*.

        If *lattice* is ``True``, the information about periodicity is printed to the resulting geo-file with ``CRYSTX`` key. If the supplied *molecule* does not contain lattice vectors (or contains less then 3 of them), this method will add them (and hence alter *molecule*!).  The length of added vectors is defined by ``default_cell_size`` class attribute.

        *settings* can also be a single string with a path to a file -- in that case this file is copied as *filename* and all the rest of this method is skipped.

        .. note::

            If *lattice* is ``True`` and the lattice present in *molecule* does not follow ReaxFF convention (the third vector aligned with Z axis, the second one with YZ plane), this method will rotate the *molecule* to fulfill these requirements.
        """
        if isinstance(settings, str) and os.path.isfile(settings):
            shutil.copy(settings, opj(self.path, filename))
        else:
            header = ['BIOGRF 200\n']

            descrp = settings.find_case('descrp')
            if descrp not in settings:
                settings.descrp = description

            for key in settings:
                val = settings[key]
                if isinstance(val, list):
                    for el in val:
                        header.append('{:6} {}\n'.format(key.upper(), el))
                elif isinstance(val, tuple):
                    header.append(('{:6} '+'{} '*len(val)+'\n').format(key.upper(), *val))
                else:
                    header.append('{:6} {}\n'.format(key.upper(), val))

            if lattice is True:
                if molecule.align_lattice(convention='reax'):
                    log("The lattice supplied for job {} did not follow the convention required by ReaxFF. I rotated the whole system for you. You're welcome".format(self._full_name()), 3)

                f = lambda x: tuple([self.default_cell_size * int(i==x) for i in range(3)])
                while len(molecule.lattice) < 3:
                    molecule.lattice.append(f(len(molecule.lattice)))

                header.append('CRYSTX  {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n'.format(*self._convert_lattice(molecule.lattice)))

            atoms = []
            for i,at in enumerate(molecule):
                newline = 'HETATM {:>5d} {:<2}               {: >10.5f}{: >10.5f}{: >10.5f} {:<2}     1 1  0.0\n'.format(i+1, at.symbol, *at.coords, at.symbol)
                atoms.append(newline)
            atoms.append('END\n')

            with open(opj(self.path, filename), 'w') as f:
                f.writelines(header)
                f.writelines(atoms)


    @staticmethod
    def _convert_lattice(lattice):
        """Convert a *lattice* expressed as three 3-dimensional vectors to (*a*, *b*, *c*, *alpha*, *beta*, *gamma*) format. Lengths of lattice vectors are expressed as *a*, *b* and *c*, angles between them as *alpha*, *beta*, *gamma*.
        """
        a, b, c = map(np.linalg.norm, lattice)
        al = np.dot(lattice[1], lattice[2])/(b*c)
        be = np.dot(lattice[0], lattice[2])/(a*c)
        ga = np.dot(lattice[0], lattice[1])/(a*b)
        al,be,ga = map(lambda x: Units.convert(np.arccos(x),'rad','deg'), [al,be,ga])
        return a, b, c, al, be, ga


#===========================================================================


def load_reaxff_control(filename, keep_order=True):
    """Return a |Settings| instance containing all data from an existing ``control`` file, indicated by *filename*.

    If *keep_order* is ``True``, the returned |Settings| instance is enriched with the ``_order`` key containing a list of all keys in the same order they were present in the loaded ``contol`` file.
    """

    if not os.path.isfile(filename):
        raise FileError('File {} not present'.format(filename))

    with open(filename, 'r') as f:
        lines = f.readlines()

    ret = Settings()
    if keep_order:
        ret._order = []

    for line in lines:
        if line.lstrip().startswith('#'):
            continue
        tmp = line.split()
        if len(tmp) > 1:
            value, key = tmp[:2]
            if len(key) > 6:
                log("load_reaxff_control: I am confused by this line:\n{}\nAre you sure {} is correct?".format(line, key), 3)
            if not value.replace('.','', 1).isdigit():
                 log("load_reaxff_control: I am confused by this line:\n{}\nAre you sure {} is correct?".format(line, value), 3)
            try:
                value = float(value) if '.' in value else int(value)
            except:
                pass
            ret[key] = value
            if keep_order:
               ret._order.append(key)
    return ret


def reaxff_control_to_settings(fpath:str) -> Settings:
    with open(fpath) as f:
        lines = [i.strip().split() for i in f]
    d = {}
    for i in lines:
        if len(i) > 1:
            try: v = int(i[0])
            except:
                try:    v = float(i[0])
                except: continue
            d[i[1]] = v

    s = Settings()
    if 'icheck' in d and d['icheck'] > 0:
        s.input.ams.task = 'singlepoint'
        if d['icheck'] == 1:
            s.input.ams.properties.gradients = True
    elif 'imetho' in d:
        v = d['imetho']
        if v == 0 or v == 2: # TODO: 2 is MD minimization, how should this be handled?
            s.input.ams.task = 'moleculardynamics'
        elif v == 1:
            s.input.ams.task = 'geometryoptimization'
        else:
            raise Exception(f"Could not understand imetho = {v}")
    else:
        raise Exception(f"Could not find a task. Is this a control file?")


    if s.input.ams.task == 'moleculardynamics':
        if 'nmdit'  in d: s.input.ams.moleculardynamics.nsteps   = d['nmdit']
        if 'tstep'  in d: s.input.ams.moleculardynamics.timestep = d['tstep']
        if 'mdtemp' in d: s.input.ams.moleculardynamics.initialvelocities.temperature = d['mdtemp']
        if 'iout2'  in d: s.input.ams.moleculardynamics.samplingfreq = d['iout2']
        if 'imdmet' in d:
            if d['imdmet'] == 1 or d['imdmet'] == 4: # NVT or NPT
                # T:
                s.input.ams.moleculardynamics.thermostat.type = 'Berendsen'
                s.input.ams.moleculardynamics.thermostat.temperature = d['mdtemp']
                if 'tdamp1' in d:
                    s.input.ams.moleculardynamics.thermostat.tau = d['tdamp1']
                if 'itdmet' in d:
                    s.input.ams.moleculardynamics.thermostat.berendsenapply = 'Local' if d['itdmet'] == 0 else 'Global'
                # P:
                if d['imdmet'] == 4: # NPT
                    s.input.ams.moleculardynamics.barostat = 'Berendsen'
                    if 'mdpres' in d:
                        s.input.ams.moleculardynamics.barostat.pressure = d['mdpres']*1e9
                    if 'pdamp1' in d:
                        s.input.ams.moleculardynamics.barostat.tau = d['pdamp1']

    elif s.input.ams.task == 'geometryoptimization':
        grads_set, fire = False, True
        if 'endmm' in d:
            s.input.ams.geometryoptimization.convergence.gradients = d['endmm'] * Units.conversion_ratio('kcal/mol','hartree')
            grads_set = True
        if 'imaxmo':
            v = d['imaxmo']
            s.input.ams.geometryoptimization.pretendconverged = True
            if v == 0:
                s.input.ams.geometryoptimization.method = 'ConjugateGradients'
                fire = False
            else:
                s.input.ams.geometryoptimization.convergence.step = v*1e-6
        if 'imaxit' in d:
            s.input.ams.geometryoptimization.pretendconverged = True
            s.input.ams.geometryoptimization.maxiterations = d['imaxit']
        if 'icelop' in d and d['icelop'] > 0:
            s.input.ams.geometryoptimization.optimizelattice = True
            fire = True # CG does not support lattice optimization

        if grads_set and fire:
            # FIRE calculates rms on atomic gradients (1/Natoms) whereas others calculate components (1/3Natoms)
            s.input.ams.geometryoptimization.convergence.gradients *= np.sqrt(3)
            s.input.ams.geometryoptimization.method = 'FIRE'

    return s
