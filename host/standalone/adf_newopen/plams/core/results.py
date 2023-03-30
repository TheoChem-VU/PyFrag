import copy
import functools
import glob
import inspect
import operator
import os
import shutil
import threading

from os.path import join as opj
from subprocess import PIPE

from .private import saferun
from .errors import ResultsError, FileError
from .functions import config, log


__all__ = ['Results']



def _caller_name_and_arg(frame):
    """Extract information about name and arguments of a function call from a *frame* object"""
    if frame is None:
        return None, None
    caller_name = frame.f_code.co_name
    caller_varnames = frame.f_code.co_varnames
    caller_arg = None
    if len(caller_varnames) > 0:
        try:
            loc = frame.f_locals
        except:
            loc = {}
        if caller_varnames[0] in loc:
            caller_arg = loc[caller_varnames[0]]
    return caller_name, caller_arg


def _privileged_access():
    """Analyze contents of the current stack to find out if privileged access to the |Results| methods should be granted.

    Privileged access is granted to two |Job| methods: |postrun| and :meth:`~scm.plams.core.basejob.Job.check`, but only if they are called from :meth:`~scm.plams.core.basejob.Job._finalize` of the same |Job| instance.
    """
    from .basejob import Job
    for frame in inspect.getouterframes(inspect.currentframe()):
        cal, arg = _caller_name_and_arg(frame[0])
        prev_cal, prev_arg = _caller_name_and_arg(frame[0].f_back)
        if cal in ['postrun', 'check'] and prev_cal == '_finalize' and arg == prev_arg and isinstance(arg, Job):
                return True
    return False


def _restrict(func):
    """Decorator that wraps methods of |Results| instances.

    Whenever decorated method is called, the status of associated job is checked. Depending of its value access to the method is granted, refused or the calling thread is forced to wait for the right :ref:`event<event-objects>` to be set.
    """

    @functools.wraps(func)
    def guardian(self, *args, **kwargs):
        if not self.job:
            raise ResultsError('Using Results not associated with any Job')

        if self.job.status in ['successful', 'copied']:
            return func(self, *args, **kwargs)

        elif self.job.status in ['preview']:
            if config.ignore_failure:
                log('WARNING: Trying to obtain results of job {} run in a preview mode. Returned value is None'.format(self.job.name), 3)
                return None
            else:
                raise ResultsError('Using Results associated with job run in a preview mode')

        elif self.job.status in ['deleted']:
            raise ResultsError('Using Results associated with deleted job')

        elif self.job.status in ['crashed', 'failed']:
            if func.__name__ == 'wait': #waiting for crashed of failed job should not trigger any warnings/exceptions
                cal, arg = _caller_name_and_arg(inspect.currentframe())
                if isinstance(arg, Results):
                    return func(self, *args, **kwargs)
            if config.ignore_failure:
                log('WARNING: Trying to obtain results of crashed or failed job {}'.format(self.job.name), 3)
                try:
                    ret = func(self, *args, **kwargs)
                except:
                    log('Obtaining results of {} failed. Returned value is None'.format(self.job.name), 3)
                    return None
                log('Obtaining results of {} successful. However, no guarantee that they make sense'.format(self.job.name), 3)
                return ret
            else:
                raise ResultsError('Using Results associated with crashed or failed job')

        elif self.job.status in ['created', 'started', 'registered', 'running']:
            log('Waiting for job {} to finish'.format(self.job.name), 1)
            if _privileged_access():
                self.finished.wait()
            else:
                self.done.wait()
            return func(self, *args, **kwargs)

        elif self.job.status in ['finished']:
            if _privileged_access():
                return func(self, *args, **kwargs)
            log('Waiting for job {} to finish'.format(self.job.name), 1)
            self.done.wait()
            return func(self, *args, **kwargs)

    return guardian



#===========================================================================
#===========================================================================
#===========================================================================



class _MetaResults(type):
    """Metaclass for |Results|. During new |Results| instance creation it wraps all methods with :func:`_restrict` decorator ensuring proper synchronization and thread safety. Methods listed in ``_dont_restrict`` as well as "magic methods" are not wrapped."""
    _dont_restrict = ['refresh', 'collect', '_clean']
    def __new__(meta, name, bases, dct):
        for attr in dct:
            if not (attr.endswith('__') and attr.startswith('__')) and callable(dct[attr]) and (attr not in _MetaResults._dont_restrict):
                dct[attr] = _restrict(dct[attr])
        return type.__new__(meta, name, bases, dct)



#===========================================================================
#===========================================================================
#===========================================================================



class Results(metaclass=_MetaResults):
    """General concrete class for job results.

    ``job`` attribute stores a reference to associated job. ``files`` attribute is a list with contents of the job folder. ``_rename_map`` is a class attribute with the dictionary storing the default renaming scheme.

    Bracket notation (``myresults[filename]``) can be used to obtain full absolute paths to files in the job folder.

    Instance methods are automatically wrapped with the "access guardian" that ensures thread safety (see |parallel|).
    """
    _rename_map = {}

    def __init__(self, job):
        self.job = job
        self.files = []
        self.finished = threading.Event()
        self.done = threading.Event()


    def refresh(self):
        """Refresh the contents of the ``files`` list. Traverse the job folder (and all its subfolders) and collect relative paths to all files found there, except files with ``.dill`` extension.

        This is a cheap and fast method that should be used every time there is a risk the contents of the job folder changed and ``files`` is no longer up-to-date. For proper working of various PLAMS elements it is crucial that ``files`` always contains up-to-date information about the contents of the job folder.

        All functions and methods defined in PLAMS that could change the state of the job folder refresh the ``files`` list, so there is no need to manually call :meth:`~Results.refresh` after, for example, :meth:`~Results.rename`. If you are implementing a new method of that kind, please don't forget about refreshing.
        """
        self.files = []
        for pth, dirs, files in os.walk(self.job.path):
            relpath = os.path.relpath(pth, self.job.path)
            self.files += [opj(relpath, x) if relpath != '.' else x for x in files]
        self.files = [x for x in self.files if not x.endswith('.dill')]


    def collect(self):
        """Collect the files present in the job folder after execution of the job is finished. This method is simply :meth:`~Results.refresh` followed by renaming according to the ``_rename_map``.

        If you wish to override this function, you have to call the parent version at the beginning.
        """
        self.refresh()
        for old, new in self.__class__._rename_map.items():
            old = old.replace('$JN', self.job.name)
            new = new.replace('$JN', self.job.name)
            if old in self.files:
                os.rename(opj(self.job.path, old), opj(self.job.path, new))
                self.files[self.files.index(old)] = new
        self.refresh()


    def wait(self):
        """wait()
        Wait for associated job to finish.

        .. technical::

            This is **not** an abstract method. It does exactly what it should: nothing. All the work is done by :func:`_restrict` decorator that is wrapped around it.
        """
        pass


    def grep_file(self, filename, pattern='', options=''):
        """grep_file(filename, pattern='', options='')
        Execute ``grep`` on a file given by *filename* and search for *pattern*.

        Additional ``grep`` flags can be passed with *options*, which should be a single string containing all flags, space separated.

        Returned value is a list of lines (strings). See ``man grep`` for details.
        """
        cmd = ['grep'] + [pattern] + options.split()
        return self._process_file(filename, cmd)


    def grep_output(self, pattern='', options=''):
        """grep_output(pattern='', options='')
        Shortcut for :meth:`~Results.grep_file` on the output file."""
        try:
            output = self.job._filename('out')
        except AttributeError:
            raise ResultsError('Job {} does not seem to be an instance of SingleJob, it does not have _filenames dictionary'.format(self.job.name))
        except KeyError:
            raise ResultsError('Job {} does not have an output'.format(self.job.name))
        return self.grep_file(output, pattern, options)


    def awk_file(self, filename, script='', progfile=None, **kwargs):
        """awk_file(filename, script='', progfile=None, **kwargs)
        Execute an AWK script on a file given by *filename*.

        The AWK script can be supplied in two ways: either by directly passing the contents of the script (should be a single string) as the *script* argument, or by providing the path (absolute or relative to *filename*) to a text file with an AWK script as the *progfile* argument. If *progfile* is not ``None``, *script* is ignored.

        Other keyword arguments (*\*\*kwargs*) can be used to pass additional variables to AWK (see ``-v`` flag in AWK manual)

        Returned value is a list of lines (strings). See ``man awk`` for details.
        """
        cmd = ['awk']
        for k,v in kwargs.items():
            cmd += ['-v', '{}={}'.format(k,v)]
        if progfile:
            if os.path.isfile(progfile):
                cmd += ['-f', progfile]
            else:
                raise FileError('File {} not present'.format(progfile))
        else:
            cmd += [script]
        return self._process_file(filename, cmd)


    def awk_output(self, script='', progfile=None, **kwargs):
        """awk_output(script='', progfile=None, **kwargs)
        Shortcut for :meth:`~Results.awk_file` on the output file."""
        try:
            output = self.job._filename('out')
        except AttributeError:
            raise ResultsError('Job {} does not seem to be an instance of SingleJob, it does not have _filenames dictionary'.format(self.job.name))
        except KeyError:
            raise ResultsError('Job {} does not have an output'.format(self.job.name))
        return self.awk_file(output, script, progfile, **kwargs)


    def rename(self, old, new):
        """rename(old, new)
        Rename a file from ``files``. In both *old* and *new* the shortcut ``$JN`` for job name can be used."""
        old = old.replace('$JN', self.job.name)
        new = new.replace('$JN', self.job.name)
        self.refresh()
        if old in self.files:
            os.rename(opj(self.job.path, old), opj(self.job.path, new))
            self.files[self.files.index(old)] = new
        else:
            raise FileError('File {} not present in {}'.format(old, self.job.path))


    def get_file_chunk(self, filename, begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None):
        """get_file_chunk(filename, begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None)

        Extract a chunk of a text file given by *filename*, consisting of all the lines between a line containing *begin* and a line containing *end*.

        *begin* and *end* should be simple strings (no regular expressions allowed) or ``None`` (in that case matching is done from the beginning or until the end of the file). If multiple blocks delimited by *begin* end *end* are present in the file, *match* can be used to indicate which one should be printed (*match*=0 prints all of them). *inc_begin* and *inc_end* can be used to include the delimiting lines in the final result (by default they are excluded).

        The returned value is a list of strings. *process* can be used to provide a function executed on each element of this list before returning it.
        """
        current_match = 0
        ret = []
        switch = (begin == None)

        append = lambda x: ret.append(x.rstrip('\n')) if (match in [0,current_match]) else None

        with open(self[filename], 'r') as f:
            for line in f:
                if switch and end and (end in line):
                    switch = False
                    if inc_end: append(line)
                    if match == current_match: break
                if switch:
                    append(line)
                if (not switch) and begin and (begin in line):
                    switch = True
                    current_match += 1
                    if inc_begin: append(line)

        return list(map(process, ret)) if process else ret


    def get_output_chunk(self, begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None):
        """get_output_chunk(begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None)
        Shortcut for :meth:`~Results.get_file_chunk` on the output file."""
        try:
            output = self.job._filename('out')
        except AttributeError:
            raise ResultsError('Job {} is not an instance of SingleJob, it does not have an output'.format(self.job.name))
        return self.get_file_chunk(output, begin, end, match, inc_begin, inc_end, process)


    def recreate_molecule(self):
        """Recreate the input molecule for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        The definiton here serves as a deafult fall-back template preventing |load_external| from crashing when a particular |Results| subclass does not define it's own :meth:`recreate_molecule`.
        """
        return None


    def recreate_settings(self):
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        The definiton here serves as a deafult fall-back template preventing |load_external| from crashing when a particular |Results| subclass does not define it's own :meth:`recreate_settings`.
        """
        return None



#=======================================================================


    def _clean(self, arg):
        """Clean the job folder. *arg* should be a string or a list of strings. See |cleaning| for details."""
        if arg == 'all':
            return

        path = self.job.path
        absfiles = [opj(path, f) for f in self.files]
        childnames = [child.name for child in self.job] if hasattr(self.job, 'children') else []
        if arg in ['none', [], None]:
            [os.remove(f) for f in absfiles if os.path.isfile(f)]

        elif isinstance(arg, list):
            rev = False
            if arg[0] == '-':
                rev = True
                arg = arg[1:]

            absarg = []
            for i in arg:
                s = i.replace('$JN', self.job.name)
                if s.find('$CH') != -1:
                    absarg += [opj(path, s.replace('$CH', ch)) for ch in childnames]
                else:
                    absarg.append(opj(path, s))

            if absarg:
                absarg = functools.reduce(operator.iadd, map(glob.glob, absarg))

            for f in absfiles:
                if (f in absarg) == rev and os.path.isfile(f):
                    os.remove(f)
                    log('Deleting file '+f, 5)

        else:
            log('WARNING: {} is not a valid keep/save argument'.format(arg), 3)
        self.refresh()


    def _copy_to(self, newresults):
        """_copy_to(newresults)
        Copy these results to *newresults*.

        This method is used when |RPM| discovers an attempt to run a job identical to a previously run job. Instead of the execution, results of the previous job are copied/linked to the new one.

        This method is called from |Results| of the old job and *newresults* should be |Results| of the new job. The goal is to faithfully recreate the state of this |Results| instance in ``newresults``. To achieve that, all the contents of the job folder are copied (or hardlinked, if your platform allows that and ``job.settings.link_files`` is ``True``) to other's job folder. Moreover, all attributes of this |Results| instance (other than ``job`` and ``files``) are exported to *newresults* using :meth:`~Results._export_attribute` method.
        """
        for name in self.files:
            newname = Results._replace_job_name(name, self.job.name, newresults.job.name)
            oldpath = opj(self.job.path, name)
            newpath = opj(newresults.job.path, newname)
            os.makedirs(os.path.dirname(newpath), exist_ok=True)
            if os.name == 'posix' and self.job.settings.link_files is True:
                os.link(oldpath, newpath)
            else:
                shutil.copy(oldpath, newpath)
            newresults.files.append(newname)
        for k,v in self.__dict__.items():
            if k in ['job', 'files', 'done', 'finished']: continue
            newresults.__dict__[k] = self._export_attribute(v, newresults)


    def _export_attribute(self, attr, other):
        """_export_attribute(attr, other)
        Export this instance's attribute to *other*. This method should be overridden in your |Results| subclass if it has some attributes that are not properly handled by :func:`python3:copy.deepcopy`.

        *other* is the |Results| instance, *attr* is the **value** of the attribute to be copied. See :meth:`SCMJob._export_attribute<scm.plams.interfaces.adfsuite.scmjob.SCMResults._export_attribute>` for an example implementation.
        """
        return copy.deepcopy(attr)


    @staticmethod
    def _replace_job_name(string, oldname, newname):
        """If *string* starts with *oldname*, maybe followed by some extension, replace *oldname* with *newname*."""
        return string.replace(oldname, newname) if (os.path.splitext(string)[0] == oldname) else string


#=======================================================================


    def __getitem__(self, name):
        """Magic method to enable bracket notation. Elements from ``files`` can be used to get absolute paths."""
        name = name.replace('$JN', self.job.name)
        if name in self.files:
            return opj(self.job.path, name)
        else:
            raise FileError('File {} not present in {}'.format(name, self.job.path))



    def _process_file(self, filename, command):
        """_process_file(filename, command)
        Skeleton for all file processing methods. Execute *command* (should be a list of strings) on *filename* and return output as a list of lines.
        """
        filename = filename.replace('$JN', self.job.name)
        if filename in self.files:
            process = saferun(command + [filename], cwd=self.job.path, stdout=PIPE)
            if process.returncode != 0:
                return []
            ret = process.stdout.decode().splitlines()
            return ret
        else:
            raise FileError('File {} not present in {}'.format(filename, self.job.path))


