import os
import threading
try:
    import dill as pickle
except ImportError:
    import pickle

from os.path import join as opj

from .basejob import MultiJob
from .errors import PlamsError, FileError
from .functions import config, log

__all__ = ['JobManager']



class JobManager:
    """Class responsible for jobs and files management.

    Every instance has the following attributes:

    *   ``foldername`` -- the working folder name.
    *   ``workdir`` -- the absolute path to the working folder.
    *   ``logfile`` -- the absolute path to the logfile.
    *   ``input`` -- the absolute path to the copy of the input file in the working folder.
    *   ``settings`` -- a |Settings| instance for this job manager (see below).
    *   ``jobs`` -- a list of all jobs managed with this instance (in order of |run| calls).
    *   ``names`` -- a dictionary with names of jobs. For each name an integer value is stored indicating how many jobs with that basename have already been run.
    *   ``hashes`` -- a dictionary working as a hash-table for jobs.

    The *path* argument should be be a path to a directory inside which the main working folder will be created. If ``None``, the directory from where the whole script was executed is used.

    The ``foldername`` attribute is initially set to the *folder* argument. If such a folder already exists, the suffix ``.002`` is appended to *folder* and the number is increased (``.003``, ``.004``...) until a non-existsing name is found. If *folder* is ``None``, the name ``plams_workdir`` is used, followed by the same procedure to find a unique ``foldername``.

    The ``settings`` attribute is directly set to the value of *settings* argument (unlike in other classes where they are copied) and it should be a |Settings| instance with the following keys:

    *   ``hashing`` -- chosen hashing method (see |RPM|).
    *   ``counter_len`` -- length of number appended to the job name in case of a name conflict.
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
            raise PlamsError('Invalid path: {}'.format(path))

        basename = os.path.normpath(folder) if folder else 'plams_workdir'
        self.foldername = basename
        n = 2
        while os.path.exists(opj(self.path, self.foldername)):
            self.foldername = basename + '.' + str(n).zfill(3)
            n += 1

        self.workdir = opj(self.path, self.foldername)
        self.logfile = opj(self.workdir, 'logfile')
        self.input = opj(self.workdir, 'input')
        os.mkdir(self.workdir)



    def load_job(self, filename):
        """Load previously saved job from *filename*.

        *Filename* should be a path to a ``.dill`` file in some job folder. A |Job| instance stored there is loaded and returned. All attributes of this instance removed before pickling are restored. That includes ``jobmanager``, ``path`` (the absolute path to the folder containing *filename* is used) and ``default_settings`` (a list containing only ``config.job``).

        See |pickling| for details.
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
                for otherjob in job.other_jobs():
                    setstate(otherjob, opj(path, otherjob.name), job)

            job.results.refresh()
            h = job.hash()
            if h is not None:
                self.hashes[h] = job
            for key in job._dont_pickle:
                job.__dict__[key] = None

        if os.path.isfile(filename):
            filename = os.path.abspath(filename)
        else:
            raise FileError('File {} not present'.format(filename))
        path = os.path.dirname(filename)
        with open(filename, 'rb') as f:
            try:
                job = pickle.load(f)
            except Exception as e:
                log("Unpickling of {} failed. Caught the following Exception:\n{}".format(filename, e), 1)
                return None

        setstate(job, path)
        return job



    def remove_job(self, job):
        """Remove *job* from the job manager. Forget its hash."""
        if job in self.jobs:
            self.jobs.remove(job)
            job.jobmanager = None
        h = job.hash()
        if h in self.hashes and self.hashes[h] == job:
            del self.hashes[h]
        if isinstance(job, MultiJob):
            for child in job:
                self.remove_job(child)
            for otherjob in job.other_jobs():
                self.remove_job(otherjob)



    def _register_name(self, job):
        """Register the name of the *job*.

        If a job with the same name was already registered, *job* is renamed by appending consecutive integers. The number of digits in the appended number is defined by the ``counter_len`` value in ``settings``.
        """

        name = job._full_name()
        if name in self.names:
            self.names[name] += 1
            job.name += '.'+ str(self.names[name]).zfill(self.settings.counter_len)
            log('Renaming job {} to {}'.format(name, job._full_name()), 3)
        else:
            self.names[name] = 1



    def _register(self, job):
        """Register the *job*. Register job's name (rename if needed) and create the job folder."""

        log('Registering job {}'.format(job.name), 7)
        job.jobmanager = self

        self._register_name(job)

        if job.path is None:
            if job.parent:
                job.path = opj(job.parent.path, job.name)
            else:
                job.path = opj(self.workdir, job.name)
        os.mkdir(job.path)

        self.jobs.append(job)
        job.status = 'registered'
        log('Job {} registered'.format(job.name), 7)



    def _check_hash(self, job):
        """Calculate the hash of *job* and, if it is not ``None``, search previously run jobs for the same hash. If such a job is found, return it. Otherwise, return ``None``"""
        h = job.hash()
        if h is not None:
            if h in self.hashes:
                prev = self.hashes[h]
                log('Job {} previously run as {}, using old results'.format(job.name, prev.name), 1)
                return prev
            else:
                self.hashes[h] = job
        return None



    def _clean(self):
        """Clean all registered jobs according to the ``save`` parameter in their ``settings``. If ``remove_empty_directories`` is ``True``,  traverse the working directory and delete all empty subdirectories."""
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

