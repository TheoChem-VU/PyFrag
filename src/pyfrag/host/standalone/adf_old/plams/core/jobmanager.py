from __future__ import unicode_literals

import glob
import os
import shutil
import threading
try:
    import dill as pickle
except ImportError:
    import pickle

from os.path import join as opj

from .common import log
from .errors import PlamsError
from .basejob import MultiJob

__all__ = ['JobManager']

class JobManager(object):
    """Class responsible for jobs and files management.

    Every instance has the following attributes:
        *   ``folder`` -- the working folder name.
        *   ``path`` -- the absolute path to the directory with the working folder.
        *   ``workdir`` -- the absolute path to the working folder (``path/folder``).
        *   ``settings`` -- a |Settings| instance for this job manager (see below).
        *   ``jobs`` -- a list of all jobs managed with this instance (in order of |run| calls).
        *   ``names`` -- a dictionary with names of jobs. For each name an integer value is stored indicating how many jobs with that name have already been run.
        *   ``hashes`` -- a dictionary working as a hash-table for jobs.

    ``path`` and ``folder`` can be adjusted with constructor arguments *path* and *folder*. If not supplied, Python current working directory and string ``plams.`` appended with PID of the current process are used.

    ``settings`` attribute is directly set to the value of *settings* argument (unlike in other classes where they are copied) and it should be a |Settings| instance with the following keys:
        *   ``hashing`` -- chosen hashing method (see |RPM|).
        *   ``counter_len`` -- length of number appended to the job name in case of name conflict.
        *   ``remove_empty_directories`` -- if ``True``, all empty subdirectories of the working folder are removed on |finish|.
    """

    def __init__(self, settings, path=None, folder=None):

        self.settings = settings
        self.jobs = []
        self.names = {}
        self.hashes = {}

        if path is None:
            self.path = os.getcwd()
        elif os.path.isdir(path):
            self.path = os.path.abspath(path)
        else:
            raise PlamsError('Invalid path: %s'%path)

        if folder is None:
            basename = 'plams.' + str(os.getpid())
            self.foldername = basename
            i = 1
            while os.path.exists(opj(self.path, self.foldername)):
                self.foldername = basename + '_' + str(i)
                i += 1
        else:
            self.foldername = os.path.normpath(folder) #normpath removes trailing /

        self.workdir = opj(self.path, self.foldername)
        self.logfile = opj(self.workdir, self.foldername+'.log')
        self.input = opj(self.workdir, self.foldername+'.inp')
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        else:
            log('WARNING: Folder %s already exists. It is strongly advised to use a fresh folder for every run. If you experience problems check config.jobmanager.jobfolder_exists setting in plams_defaults.py' % self.workdir, 1)


    def _register_name(self, job):
        """Register the name of the *job*.

        If a job with the same name was already registered, *job* is renamed by appending consecutive integers. Number of digits in the appended number is defined by ``counter_len`` value in job manager's ``settings``.
        """

        if job.name in self.names:
            self.names[job.name] += 1
            newname = job.name +'.'+ str(self.names[job.name]).zfill(self.settings.counter_len)
            log('Renaming job %s to %s' % (job.name,newname), 3)
            job.name = newname
        else:
            self.names[job.name] = 1


    def _register(self, job):
        """Register the *job*. Register job's name (rename if needed) and create the job folder."""

        log('Registering job %s' % job.name, 7)
        job.jobmanager = self

        self._register_name(job)

        if job.path is None:
            if job.parent:
                job.path = opj(job.parent.path, job.name)
            else:
                job.path = opj(self.workdir, job.name)
        if os.path.exists(job.path):
            if self.settings.jobfolder_exists == 'remove':
                shutil.rmtree(job.path)
            elif self.settings.jobfolder_exists == 'rename':
                i = 1
                while os.path.exists(job.path + '.old' + str(i)):
                    i += 1
                newname = job.path + '.old' + str(i)
                os.rename(job.path, newname)
                log('Folder %s already present. Renaming it to %s'%(job.path, newname), 1)
            else:
                raise PlamsError('Folder %s already present in the filesystem. Consider using a fresh working folder or adjusting config.jobmanager.jobfolder_exists'%job.path)
        os.mkdir(job.path)

        self.jobs.append(job)
        job.status = 'registered'
        log('Job %s registered' % job.name, 7)


    def _check_hash(self, job):
        """Calculate the hash of *job* and, if it is not ``None``, search previously run jobs for the same hash. If such a job is found, return it. Otherwise, return ``None``"""
        h = job.hash()
        if h is not None:
            if h in self.hashes:
                prev = self.hashes[h]
                log('Job %s previously run as %s, using old results' % (job.name, prev.name), 1)
                return prev
            else:
                self.hashes[h] = job
        return None


    def load_job(self, filename):
        """Load previously saved job from *filename*.

        *Filename* should be a path to ``.dill`` file in some job folder. A |Job| instance stored there is loaded and returned. All attributes of this instance removed before pickling are restored. This includes ``jobmanager``, ``path`` (absolute path to *filename* is used), ``default_setting`` (list containing only ``config.job``) and also ``parent`` in case of children of some |MultiJob|.

        See :ref:`pickling` for details.
        """
        def setstate(job, path, parent=None):
            job.parent = parent
            job.jobmanager = self
            job.default_settings = [config.job]
            job.path = path
            if isinstance(job, MultiJob):
                job._lock = threading.Lock()
                for child in job:
                    setstate(child, opj(path, child.name), job)
            job.results.refresh()
            h = job.hash()
            if h is not None:
                self.hashes[h] = job
            for key in job._dont_pickle:
                job.__dict__[key] = None

        filename = os.path.abspath(filename)
        path = os.path.dirname(filename)
        with open(filename, 'rb') as f:
            job = pickle.load(f)

        setstate(job, path)

        return job

    def remove_job(self, job):
        """Remove *job* from job manager. Forget its hash."""
        if job in self.jobs:
            self.jobs.remove(job)
            job.jobmanager = None
        h = job.hash()
        if h in self.hashes and self.hashes[h] == job:
            del self.hashes[h]


    def _clean(self):
        """Clean all registered jobs according to their ``save`` parameter in their ``settings``. If ``remove_empty_directories`` is ``True``,  traverse the working directory and delete all empty subdirectories.
        """

        log('Cleaning job manager', 7)

        for job in self.jobs:
            job.results._clean(job.settings.save)

        if self.settings.remove_empty_directories:
            for root, dirs, files in os.walk(self.workdir, topdown=False):
                for dirname in dirs:
                    fullname = opj(root, dirname)
                    if not os.listdir(fullname):
                        os.rmdir(fullname)

        log('Job manager cleaned', 7)

