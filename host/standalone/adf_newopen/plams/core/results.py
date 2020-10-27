from __future__ import unicode_literals
from six import add_metaclass

import copy
import functools
import glob
import inspect
import operator
import os
import shutil
import threading
import time
import types
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

from os.path import join as opj

from .common import log, string
from .errors import ResultsError, FileError

__all__ = ['Results']


#===================================================================================================
#===================================================================================================
#===================================================================================================


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

    Privileged access is granted to two |Job| methods: |postrun| and :meth:`~scm.plams.basejob.Job.check`, but only if they are called from :meth:`~scm.plams.basejob.Job._finalize` of the same |Job| instance.
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
            raise ResultsError('Using Results not associated with any job')

        if self.job.status in ['successful', 'copied']:
            return func(self, *args, **kwargs)

        elif self.job.status in ['created', 'preview']:
            if config.ignore_failure:
                log("WARNING: Trying to obtain results of job %s with status '%s'. Returned value is None" % (self.job.name, self.job.status), 3)
                return None
            else:
                raise ResultsError('Using Results associated with unfinished job')

        elif self.job.status in ['crashed', 'failed']:
            if func.__name__ == 'wait': #waiting for crashed of failed job should not trigger any warnings/exceptions
                cal, arg = _caller_name_and_arg(inspect.currentframe())
                if isinstance(arg, Results):
                    return func(self, *args, **kwargs)
            if config.ignore_failure:
                log('WARNING: Trying to obtain results of crashed or failed job %s' % self.job.name, 3)
                try:
                    ret = func(self, *args, **kwargs)
                except:
                    log('Obtaining results of %s failed. Returned value is None' % self.job.name, 3)
                    return None
                log('Obtaining results of %s successful. However, no guarantee that they make sense'% self.job.name, 3)
                return ret
            else:
                raise ResultsError('Using Results associated with crashed or failed job')

        elif self.job.status in ['started', 'registered', 'running']:
            log('Waiting for job %s to finish' % self.job.name, 3)
            if _privileged_access():
                self.finished.wait()
            else:
                self.done.wait()
            return func(self, *args, **kwargs)

        elif self.job.status in ['finished']:
            if _privileged_access():
                return func(self, *args, **kwargs)
            log('Waiting for job %s to finish' % self.job.name, 3)
            self.done.wait()
            return func(self, *args, **kwargs)

    return guardian



#===================================================================================================
#===================================================================================================
#===================================================================================================


class _MetaResults(type):
    """Metaclass for |Results|. During new |Results| instance creation it wraps all methods with :func:`_restrict` decorator ensuring proper synchronization and thread safety. Methods listed in ``_dont_restrict`` as well as "magic methods" are not wrapped."""
    _dont_restrict = ['refresh', 'collect', '_clean']
    def __new__(meta, name, bases, dct):
        for attr in dct:
            if not (attr.endswith('__') and attr.startswith('__')) and callable(dct[attr]) and (attr not in _MetaResults._dont_restrict):
                dct[attr] = _restrict(dct[attr])
        return type.__new__(meta, name, bases, dct)


#===================================================================================================
#===================================================================================================
#===================================================================================================

@add_metaclass(_MetaResults)
class Results(object):
    """General concrete class for job results.

    ``job`` attribute stores a reference to associated job. ``files`` attribute is a list with contents of the job folder. ``_rename_map`` is a class attribute with the dictionary storing the default renaming scheme.

    Bracket notation (``myresults[filename]`` can be used to obtain full absolute paths to files in the job folder.

    Instance methods are automatically wrapped with access guardian which ensures thread safety (see :ref:`parallel`).
    """
    _rename_map = {}

    def __init__(self, job):
        self.job = job
        self.files = []
        self.finished = threading.Event()
        self.done = threading.Event()


    def refresh(self):
        """Refresh the contents of ``files`` list. Traverse the job folder (and all its subfolders) and collect relative paths to all files found there, except files with ``.dill`` extension.

        This is a cheap and fast method that should be used every time there is some risk that contents of the job folder changed and ``files`` list is no longer up-to-date. For proper working of various PLAMS elements it is crucial that ``files`` always contains up-to-date information about contents of job folder.

        All functions and methods defined in PLAMS that could change the state of job folder take care about refreshing ``files``, so there is no need to manually call :meth:`~Results.refresh` after, for example, :meth:`~Results.rename`. If you are implementing new method of that kind, don't forget about refreshing.
        """
        self.files = []
        for pth, dirs, files in os.walk(self.job.path):
            relpath = os.path.relpath(pth, self.job.path)
            self.files += [opj(relpath, x) if relpath != '.' else x for x in files]
        self.files = [x for x in self.files if not x.endswith('.dill')]


    def collect(self):
        """Collect the files present in the job folder after execution of the job is finished. This method is simply :meth:`~Results.refresh` plus rename according to ``_rename_map``.

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


    def awk_file(self, filename, script='', progfile=None, **kwargs):
        """awk_file(filename, script='', progfile=None, **kwargs)
        Execute an AWK script on a file given by *filename*.

        The AWK script can be supplied in two ways: either by directly passing the contents of the script (should be a single string) as a *script* argument, or by providing the path (absolute or relative to the file pointed by *filename*) to some external file containing the actual AWK script using *progfile* argument. If *progfile* is not ``None``, the *script* argument is ignored.

        Other keyword arguments (*\*\*kwargs*) can be used to pass additional variables to AWK (see ``-v`` flag in AWK manual)

        Returned value is a list of lines (strings). See ``man awk`` for details.
        """
        cmd = ['awk']
        for k,v in kwargs.items():
            cmd += ['-v', '%s=%s'%(k,v)]
        if progfile:
            if os.path.isfile(progfile):
                cmd += ['-f', progfile]
            else:
                raise FileError('File %s not present' % progfile)
        else:
            cmd += [script]
        return self._process_file(filename, cmd)


    def grep_output(self, pattern='', options=''):
        """grep_output(pattern='', options='')
        Shortcut for :meth:`~Results.grep_file` on the output file."""
        try:
            output = self.job._filename('out')
        except AttributeError:
            raise ResultsError('Job %s is not an instance of SingleJob, it does not have an output' % self.job.name)
        return self.grep_file(output, pattern, options)


    def awk_output(self, script='', progfile=None, **kwargs):
        """awk_output(script='', progfile=None, **kwargs)
        Shortcut for :meth:`~Results.awk_file` on the output file."""
        try:
            output = self.job._filename('out')
        except AttributeError:
            raise ResultsError('Job %s is not an instance of SingleJob, it does not have an output' % self.job.name)
        return self.awk_file(output, script, progfile, **kwargs)


    def rename(self, old, new):
        """rename(old, new)
        Rename a file from ``files``. In both *old* and *new* shortcut ``$JN`` can be used."""
        old = old.replace('$JN', self.job.name)
        new = new.replace('$JN', self.job.name)
        self.refresh()
        if old in self.files:
            os.rename(opj(self.job.path, old), opj(self.job.path, new))
            self.files[self.files.index(old)] = new
        else:
            raise FileError('File %s not present in %s' % (old, self.job.path))


#===============================================================================================


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
            log('WARNING: %s is not a valid keep/save argument' % str(arg), 3)
        self.refresh()


    def _copy_to(self, other):
        """_copy_to(other)
        Copy these results to *other*.

        This method is used when |RPM| discovers an attempt to run a job identical to the one previously run. Instead of execution, results of the previous job are copied/linked to the new one.

        This method is called from results of old job and *other* should be results of new job. The goal is to faithfully recreate the state of ``self`` in ``other``. To achieve that all contents of jobs folder are copied (or hardlinked, if your platform allows that and ``self.settings.link_files`` is ``True``) to other's job folder. Moreover, all attributes of ``self`` (other than ``job`` and ``files``) are exported to *other* using :meth:`~Results._export_attribute` method.
        """
        for name in self.files:
            newname = Results._replace_job_name(name, self.job.name, other.job.name)
            args = (opj(self.job.path, name), opj(other.job.path, newname))
            if os.name == 'posix' and self.job.settings.link_files is True:
                os.link(*args)
            else:
                shutil.copy(*args)
            other.files.append(newname)
        for k,v in self.__dict__.items():
            if k in ['job', 'files', 'done', 'finished']: continue
            other.__dict__[k] = self._export_attribute(v, other)


    def _export_attribute(self, attr, other):
        """_export_attribute(attr, other)
        Export this instance's attribute to *other*. This method should be overridden in your |Results| subclass if it has some attribute that is not properly copyable by :func:`python2:copy.deepcopy`.

        *other* is the |Results| instance, *attr* is the **value** of the attribute to be copied. See :meth:`SCMJob._export_attribute<scm.plams.scmjob.SCMResults._export_attribute>` for an example implementation.
        """
        return copy.deepcopy(attr)


    @staticmethod
    def _replace_job_name(string, oldname, newname):
        """If *string* starts with *oldname*, maybe followed by some extension, replace *oldname* with *newname*."""
        return string.replace(oldname, newname) if (os.path.splitext(string)[0] == oldname) else string

    #===============================================================================================

    def __getitem__(self, name):
        """Magic method to enable bracket notation. Elements from ``files`` can be used to get absolute paths."""
        name = name.replace('$JN', self.job.name)
        if name in self.files:
            return opj(self.job.path, name)
        else:
            raise FileError('File %s not present in %s' % (name, self.job.path))



    def _process_file(self, filename, command):
        """_process_file(filename, command)
        Skeleton for all file processing methods. Execute *command* (should be a list of strings) on *filename* and return output as a list of lines.
        """
        filename = filename.replace('$JN', self.job.name)
        if filename in self.files:
            try:
                output = subprocess.check_output(command + [filename], cwd=self.job.path)
            except subprocess.CalledProcessError:
                return []
            output = string(output)
            ret = output.split('\n')
            if ret[-1] == '':
                ret = ret[:-1]
            return ret
        else:
            raise FileError('File %s not present in %s' % (filename, self.job.path))


