from __future__ import unicode_literals

import os
import shutil
import struct
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

from bisect import bisect

from ..core.errors import FileError
from ..core.common import log

__all__ = ['KFFile', 'KFReader']

class KFReader(object):
    """A class for efficient Python-native reader of binary files in KF format.

    This class offers read-only access to any fragment of data from a KF file. Unlike other Python KF readers, this one does not use the Fortran binary ``dmpkf`` to process KF files, but instead reads and interprets raw binary data straight from the file, on Python level. That approach results in significant speedup (by a factor of few hundreds for large files extracted variable by variable).

    The constructor argument *path* should be a string with a path (relative or absolute) to an existing KF file.

    *blocksize* indicates the length of basic KF file block. So far, all KF files produced by any of ADFSuite programs have the same block size of 4096 bytes. Unless you're doing something *very* special, you should not touch this value.

    Organization of data inside KF file can depend on a machine on which this file was produced. Two parameters can vary: the length of integer (32 or 64 bit) and endian (little or big). These parameters have to be determined before any reading can take place, otherwise the results will have no sense. If the constructor argument *autodetect* is ``True``, the constructor attempts to automatically detect the format of a given KF file, allowing to read files created on a machine with different endian or integer length. This automatic detection is enabled by default and it is advised to leave it that way. If you wish to disable it, you should set ``endian`` and ``word`` attributes manually before reading anything (see the code for details).

    .. note ::

        This class consists of quite technical, low level code. If you don't need to modify or extend |KFReader|, you can safely ignore all private methods, all you need is :meth:`~KFReader.read` and occasionally :meth:`~KFReader.__iter__`

    """

    _sizes = {'s':1,'i':4,'d':8,'q':8}

    def __init__(self, path, blocksize=4096, autodetect=True):
        if os.path.isfile(path):
            self.path = os.path.abspath(path)
        else:
            raise FileError('No such file: %s' % path)

        self._blocksize = blocksize
        self.endian = '<'   # endian: '<' = little, '>' = big
        self.word = 'i'     # length of int: 'i' = 4 bits, 'q' = 8 bits
        self._sections = None
        if autodetect:
            self._autodetect()


    def read(self, section, variable):
        """Extract and return data for a *variable* located in a *section*.

        For single-value numerical or boolean variables returned value is a single number or bool. For longer variables this method returns a list of values. For string variables a single string is returned.
        """

        if self._sections is None:
            self._create_index()

        try:
            tmp = self._sections[section]
        except KeyError:
            log('KFReader.read: Section %s not present in %s. Returning None as value' % (section, self.path), 1)
            return None
        try:
            vtype, vlb, vstart, vlen = tmp[variable]
        except KeyError:
            log('KFReader.read: Variable %s not present in section %s of %s. Returning None as value' % (variable, section, self.path))
            return None

        ret = []
        first = True
        with open(self.path, 'rb') as f:
            for i in KFReader._datablocks(self._data[section], vlb):
                if first:
                    ret = self._get_data(self._read_block(f,i))[vtype-1][vstart-1:]
                    first = False
                else:
                    ret += self._get_data(self._read_block(f,i))[vtype-1]
                if len(ret) >= vlen:
                    ret = ret[:vlen]
                    if len(ret) == 1: return ret[0]
                    return ret


    def __iter__(self):
        """Iteration yields pairs of section name and variable name."""
        if self._sections is None:
            self._create_index()
        for section in self._sections:
            for variable in self._sections[section]:
                yield section, variable


    def _settings_reduce(self):
        """When this object is present as a value in some |Settings| instance and string representation is needed, use the absolute path. See :meth:`Settings.__reduce__<scm.plams.settings.Settings.__reduce__>` for details."""
        return self.path


    def _autodetect(self):
        """Try to automatically detect the format (int size and endian) of this KF file."""
        with open(self.path, 'rb') as f:
            b = f.read(128)

        one = b[80:84]

        if struct.unpack(b'32s',b[48:80])[0] == b'SUPERINDEX                      ':
            self.word = 'i'
        elif struct.unpack(b'32s',b[64:96])[0] == b'SUPERINDEX                      ':
            self.word = 'q'
            one = b[96:104]
        else:
            log('KFFile._autodetect: Unable to detect integer size and endian of %s. Using defaults (4 bytes and little endian)' % self.path, 3)
            return

        for e in ['<', '>']:
            if struct.unpack(str(e+self.word), one)[0] == 1:
                self.endian = e
                d = {'q':'8 bytes','i':'4 bytes','<':'little endian','>':'big endian',}
                log(('Format of {0} detected to {'+self.word+'} and {'+self.endian+'}').format(self.path, **d), 7)


    def _read_block(self, f, pos):
        """Read a single block of binary data from posistion *pos* in file *f*."""
        f.seek((pos-1)*self._blocksize)
        return f.read(self._blocksize)


    def _parse(self, block, format):   #format = [(32,'s'),(4,'i'),(2,'d')]
        """Translate a *block* of binary data into list of values in specified *format*.

        *format* should be a list of pairs *(a,t)* where *t* is one of the following characters: ``'s'`` for string, ``'i'`` for 32-bit integer, ``'q'`` for 64-bit integer and *a* is the number of occurrences (or length of a string).

        For example, if *format* is equal to ``[(32,'s'),(4,'i'),(2,'d'),(2,'i')]``, the contents of *block* are divided into 72 bytes (32*1 + 4*4 + 2*8 + 2*4 = 72) chunks (possibly droping the last one, if it's shorter than 72 bytes). Then each chunk is translated to a 9-tuple of string, 4 ints, 2 floats and 2 ints. List of such tuples is the returned value.
        """
        step = 0
        formatstring = self.endian
        for a,t in format:
            step += a * self._sizes[t]
            formatstring += str(a) + t

        ret = []
        if step > 0:
            while step <= len(block):
                new = struct.unpack(str(formatstring), block[:step])
                new = tuple(map(lambda x: x.decode() if isinstance(x,bytes) else x, new))
                ret.append(new)
                block = block[step:]
        return ret


    def _get_data(self, datablock):
        """Extract all data from a single data block. Returned value is a 4-tuple of lists, one list for each data type (respectively: int, float, str, bool).
        """
        hlen = 4 * self._sizes[self.word]
        i,d,s,b = self._parse(datablock[:hlen],[(4,self.word)])[0]
        contents = self._parse(datablock[hlen:], zip((i,d,s,b),(self.word,'d','s',self.word)))
        if contents:
            ret = list(contents[0])
            return ret[:i], ret[i:i+d], ret[i+d], map(bool,ret[i+d+1:])
        else:
            return [], [], [], []


    def _create_index(self):
        """Find and parse relevant index blocks of KFFile to extract the information about location of all sections and variables.

        Two dictionaries are populated during this process. ``_data`` contains, for each section, a list of triples describing how logical blocks of data are mapped into physical ones. For example, ``_data['General'] = [(3,6,12), (9,40,45)]`` means that logical blocks 3-8 of section ``General`` are located in physical blocks 6-11 and logical blocks 9-13 in physical blocks 40-44. This list is always sorted via first tuple elements allowing efficient access to arbitrary logical block of each section.

        The second dictionary, ``_sections``, is used to locate each variable within its section. For each section, it contains another dictionary of each variable of this section. So ``_section[sec][var]`` contains all information needed to extract variable ``var`` from section ``sec``. This is a 4-tuple containing the following information: variable type, logic block in which the variable first occurs, position within this block where its data start and the length of the variable. Combining this information with mapping stored in ``_data`` allows to extract each single variable.
        """

        hlen = 32 + 7 * self._sizes[self.word]   #length of index block header

        with open(self.path, 'rb') as f:
            superlist = self._parse(self._read_block(f, 1), [(32,'s'),(4,self.word)])
            nextsuper = superlist[0][4]
            while nextsuper != 1:
                nsl = self._parse(self._read_block(f, nextsuper), [(32,'s'),(4,self.word)])
                nextsuper = nsl[0][4]
                superlist += nsl

            self._data = {}   #list of triples to convert logical to physical block numbers
            self._sections = {}
            for key, pb, lb, le, ty in superlist:   #pb=physical block, lb=logical block, le=length, ty=type (3 for index, 4 for data)
                key = key.rstrip(' ')
                if key in ['SUPERINDEX', 'EMPTY']: continue
                if ty == 4:   #data block
                    if key not in self._data:
                        self._data[key] = []
                    self._data[key].append((lb, pb, pb+le))
                elif ty == 3:   #index block
                    if key not in self._sections:
                        self._sections[key] = {}
                    for i in range(le):
                        indexblock = self._read_block(f, pb+i)
                        header = self._parse(indexblock[:hlen],[(32,'s'),(7,self.word)])[0]
                        body = self._parse(indexblock[hlen:],[(32,'s'),(6,self.word)])
                        for var, vlb, vstart, vlen, _xx1, _xx2, vtype in body:
                            var = var.rstrip(' ')
                            if var == 'EMPTY': continue
                            self._sections[key][var] = (vtype, vlb, vstart, vlen)

            for k,v in self._data.items():
                self._data[k] = sorted(v)


    @staticmethod
    def _datablocks(lst, n=1):
        """Transform a list of tuples ``[(x1,a1,b1),(x2,a2,b2),...]`` into an iterator over ``range(a1,b1)+range(a2,b2)+...`` Iteration starts from nth element of this list."""
        i = bisect(list(zip(*lst))[0], n) - 1
        lb, first, last = lst[i]
        ret = first + n - lb
        while i < len(lst):
            while ret < last:
                yield ret
                ret += 1
            i += 1
            if i < len(lst):
                _, ret, last = lst[i]


#===================================================================================================
#===================================================================================================

class KFFile(object):
    """A class for reading and writing binary files in KF format.

    This class acts as a wrapper around |KFReader| collecting all the data written by user in some "temporary zone" and using Fortran binaries ``udmpkf`` and ``cpkf`` to write this data to the physical file when needed.

    The constructor argument *path* should be a string with a path to an existing KF file or a new KF file that you wish to create. If a path to existing file is passed, new |KFReader| instance is created allowing to read all the data from this file.

    When :meth:`~KFFile.write` method is used, the new data is not immediately written to a disk. Instead of that, it is temporarily stored in ``tmpdata`` dictionary. When method :meth:`~KFFile.save` is invoked, contents of that dictionary are written to a physical file and ``tmpdata`` is emptied.

    Other methods like :meth:`~KFFile.read` or :meth:`~KFFile.delete_section` are aware of ``tmpdata`` and work flawlessly, regardless if :meth:`~KFFile.save` was called or not.

    By default, :meth:`~KFFile.save` is automatically invoked after each :meth:`~KFFile.write`, so physical file on a disk is always "actual". This behavior can be adjusted with *autosave* constructor parameter. Having autosave enabled is usually a good idea, however, if you need to write a lot of small pieces of data to your file, the overhead of calling ``udmpkf`` and ``cpkf`` after *every* :meth:`~KFFile.write` can lead to significant delays. In such a case it is advised to disable autosave and call :meth:`~KFFile.save` manually, when needed.

    Dictionary-like bracket notation can be used as a shortcut to read and write variables::

        mykf = KFFile('someexistingkffile.kf')
        #all three below are equivalent
        x = mykf['General%Termination Status']
        x = mykf[('General','Termination Status')]
        x = mykf.read('General','Termination Status')

        #all three below are equivalent
        mykf['Geometry%xyz'] = somevariable
        mykf[('Geometry','xyz')] = somevariable
        mykf.write('Geometry','xyz', somevariable)

    """
    _types = {int : (1, 8, lambda x:'%10i'%x),
            float : (2, 3, lambda x:'%26.16e'%x),
              str : (3, 80, lambda x: x),
             bool : (4, 80, lambda x: 'T' if x else 'F')}


    def __init__(self, path, autosave=True):
        self.autosave = autosave
        self.path = os.path.abspath(path)
        self.tmpdata = {}
        self.reader = KFReader(self.path) if os.path.isfile(self.path) else None


    def read(self, section, variable):
        """Extract and return data for a *variable* located in a *section*.

        For single-value numerical or boolean variables returned value is a single number or bool. For longer variables this method returns a list of values. For string variables a single string is returned.
        """
        if section in self.tmpdata and variable in self.tmpdata[section]:
            return self.tmpdata[section][variable]
        return self.reader.read(section, variable)


    def write(self, section, variable, value):
        """Write a *variable* with a *value* in a *section* . If such a variable already exists in this section, the old value is overwritten."""
        if not isinstance(value, (int,bool,float,str,list)):
            raise ValueError('Trying to store improper value in KFFile')
        if isinstance(value, list):
            if len(value) == 0:
                raise ValueError('Cannot store empty lists in KFFile')
            if any(not isinstance(i, type(value[0])) for i in value):
                raise ValueError('Lists stored in KFFile must have all elements of the same type')
            if not isinstance(value[0], (int,bool,float)):
                raise ValueError('Only lists of int, float or bool can be stored in KFFile')

        if section not in self.tmpdata:
            self.tmpdata[section] = {}
        self.tmpdata[section][variable] = value

        if self.autosave:
            self.save()


    def save(self):
        """Save all changes stored in ``tmpdata`` to physical file on a disk."""
        if len(self.tmpdata) > 0 and any(len(i) > 0 for i in self.tmpdata.values()):
            txt = ''
            newvars = []
            for section in self.tmpdata:
                for variable in self.tmpdata[section]:
                    val = self.tmpdata[section][variable]
                    txt += '%s\n%s\n%s\n' % (section, variable, KFFile._str(val))
                    newvars.append(section+'%'+variable)
            self.tmpdata = {}

            tmpfile = self.path+'.tmp' if self.reader else self.path
            with open(os.devnull, 'wb') as null:
                popen = subprocess.Popen(['udmpkf', tmpfile], stdin=subprocess.PIPE, stdout=null, stderr=null)
                popen.communicate(txt.encode())

            if self.reader:
                with open(os.devnull, 'wb') as null:
                    subprocess.call(['cpkf', tmpfile, self.path] + newvars, stdout=null, stderr=null)
                os.remove(tmpfile)
            self.reader = KFReader(self.path)


    def delete_section(self, section):
        """Delete the entire *section* from this KF file."""
        if section in self.tmpdata:
            del self.tmpdata[section]
        if self.reader:
            if not self.reader._sections:
                self.reader._create_index()
            if section in self.reader._sections:
                tmpfile = self.path+'tmp'
                with open(os.devnull, 'wb') as null:
                    subprocess.call(['cpkf', self.path, tmpfile, '-rm', section], stdout=null, stderr=null)
                shutil.move(tmpfile, self.path)
                self.reader = KFReader(self.path)



    def __getitem__(self, name):
        """Allow to use ``x = mykf['section%variable']`` or ``x = mykf[('section','variable')]`` instead of ``x = kf.read('section', 'variable')``."""
        section, variable = KFFile._split(name)
        return self.read(section, variable)

    def __setitem__(self, name, value):
        """Allow to use ``mykf['section%variable'] = value`` or ``mykf[('section','variable')] = value`` instead of ``kf.write('section', 'variable', value)``."""
        section, variable = KFFile._split(name)
        self.write(section, variable, value)

    def __iter__(self):
        """Iteration yields pairs of section name and variable name."""
        ret = []
        if self.reader:
            for sec,var in self.reader:
                ret.append((sec,var))
        for sec in self.tmpdata:
            for var in self.tmpdata[sec]:
                ret.append((sec,var))
        ret.sort(key=lambda x: x[0]+x[1])
        for i in ret:
            yield i



    def _settings_reduce(self):
        """When this object is present as a value in some |Settings| instance and string representation is needed, use the absolute path. See :meth:`Settings.__reduce__<scm.plams.settings.Settings.__reduce__>` for details."""
        return self.path



    @staticmethod
    def _split(name):
        """Ensure that a key used in bracket notation is of the form ``'section%variable'`` or ``('section','variable')``. If so, return a tuple ``('section','variable')``."""
        if isinstance(name, tuple) and len(name) == 2 and isinstance(name[0], str) and isinstance(name[1], str):
                return name[0], name[1]
        if isinstance(name, str):
            s = name.split('%')
            if len(s) == 2:
                return s[0], s[1]
        raise ValueError('Improper key used in KFFile dictionary-like notation')


    @staticmethod
    def _str(val):
        """Return a string representation of *val* in the form that can be understood by ``udmpkf``."""
        if isinstance(val, (int,bool,float)):
            val = [val]
        t,step,f = KFFile._types[type(val[0])]
        l = len(val) #length for str == 0 mod 160?
        ret = '%10i%10i%10i'%(l,l,t)
        for i,el in enumerate(val):
            if i%step == 0: ret += '\n'
            ret += f(el)
        return ret
