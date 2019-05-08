from __future__ import unicode_literals

import copy
import hashlib
import os
import stat
import threading
import time
import types
try:
    import dill as pickle
except ImportError:
    import pickle

from os.path import join as opj

from .common import log
from .errors import PlamsError, ResultsError
from .results import Results
from .settings import Settings

__all__ = ['SingleJob','MultiJob']


class Job(object):
    """General abstract class for all kind of computational tasks.

    Methods common for all kinds of jobs are gathered here. Instances of |Job| should never be created. It should not be subclassed either. If you wish to define a new type of job please subclass either |SingleJob| or |MultiJob|.

    Methods that are meant to be explicitly called by the user are |run| and occasionally :meth:`~Job.pickle`. In most cases |pickling| is done automatically, but if for some reason you wish to do it manually, you can use :meth:`~Job.pickle` method.

    Methods that can be safely overridden in subclasses are:
        *   :meth:`~Job.check`
        *   :meth:`~Job.hash` (see |RPM|)
        *   |prerun| and |postrun| (see :ref:`prerun-postrun`)

    Other methods should remain unchanged.

    Class attribute ``_result_type`` defines the type of results associated with this job. It should point to a class and it **must** be a |Results| subclass.

    Every job instance has the following attributes. Values of these attributes are adjusted automatically and should not be set by the user:
        *   ``status`` -- current status of the job in human-readable format.
        *   ``results`` -- reference to a results instance. An empty instance of the type stored in ``_result_type`` is created when the job constructor is called.
        *   ``path`` -- an absolute path to the job folder.
        *   ``jobmanager`` -- a job manager associated with this job.
        *   ``parent`` -- a pointer to the parent job if this job is a child job of some |MultiJob|. ``None`` otherwise.

    These attributes can be modified, but only before |run| is called:
        *   ``name`` -- the name of the job.
        *   ``settings`` -- settings of the job.
        *   ``default_settings`` -- see :ref:`default-settings`.
        *   ``depend`` -- a list of explicit dependencies.
        *   ``_dont_pickle`` -- additional list of this instance's attributes that will be removed before pickling. See :ref:`pickling` for details.

    """

    _result_type = Results

    def __init__(self, name='plamsjob', settings=None, depend=None):
        if os.path.sep in name:
            raise PlamsError('Job name cannot contain %s'%os.path.sep)
        self.status = 'created'
        self.results = self.__class__._result_type(self)
        self.name = name
        self.path = None
        self.jobmanager = None
        self.parent = None
        self.settings = Settings()
        self.default_settings = [config.job]
        self.depend = depend or []
        self._dont_pickle = []
        if settings is not None:
            if isinstance(settings, Settings):
                self.settings = settings.copy()
            if isinstance(settings, Job):
                self.settings = settings.settings.copy()



    def __getstate__(self):
        """Prepare an instance for pickling.

        Attributes ``jobmanager``, ``parent``, ``default_settings`` and ``_lock`` are removed, as well as all attributes listed in ``self._dont_pickle``.
        """
        remove = ['jobmanager', 'parent', 'default_settings', '_lock'] + self._dont_pickle
        return {k:v for k,v in self.__dict__.items() if k not in remove}



    def run(self, jobrunner=None, jobmanager=None, **kwargs):
        """Run the job using *jobmanager* and *jobrunner* (or defaults, if ``None``). Other keyword arguments (*\*\*kwargs*) are stored in ``run`` branch of job's settings. Returned value is the |Results| instance associated with this job.

        .. note::

            This method should **not** be overridden.

        .. technical::

            This method does not do too much by itself. After simple initial preparation it passes control to job runner, which decides if a new thread should be started for this job. The role of the job runner is to execute three methods that make the full job life cycle: :meth:`~Job._prepare`, :meth:`~Job._execute` and :meth:`~Job._finalize`. During :meth:`~Job._execute` the job runner is called once again to execute the runscript (only in case of |SingleJob|).
        """
        if self.status != 'created':
            raise PlamsError('Trying to run previously started job %s' % self.name)

        self.status = 'started'
        log('Job %s started' % self.name, 1)

        self.settings.run.soft_update(Settings(kwargs))

        jobrunner = jobrunner or config.default_jobrunner
        jobmanager = jobmanager or config.jm

        jobrunner._run_job(self, jobmanager)
        return self.results



    def pickle(self, filename=None):
        """Pickle this instance and save to a file indicated by *filename*. If ``None``, save to ``[jobname].dill`` in the job folder."""
        filename = filename or opj(self.path, self.name+'.dill')
        with open(filename, 'wb') as f:
            try:
                pickle.dump(self, f, -1)
            except:
                log("Pickling of %s failed" % self.name, 1)



    def check(self):
        """Check if the calculation was successful.

        This method can be overridden in concrete subclasses for different types of jobs. It should return a boolean value.

        The definition here serves as a default, to prevent crashing if a subclass does not define its own :meth:`~scm.plams.basejob.Job.check`. It always returns ``True``.
        """
        return True



    def hash(self):
        """Calculate the hash of this instance. Abstract method."""
        raise PlamsError('Trying to run an abstract method Job.hash()')



    def prerun(self):
        """Actions to take before the actual job execution.

        This method is initially empty, it can be defined in subclasses or directly added to either whole class or a single instance using |binding_decorators|.
        """
        pass

    def postrun(self):
        """Actions to take just after the actual job execution.

        This method is initially empty, it can be defined in subclasses or directly added to either whole class or a single instance using |binding_decorators|.
        """
        pass

    #===============================================================================================


    def _prepare(self, jobmanager):
        """Prepare the job for execution. This method collects steps 1-7 from :ref:`job-life-cycle`. Should not be overridden. Returned value indicates if job execution should continue (|RPM| did not find this job previously run)."""

        log('Starting %s._prepare()' % self.name, 7)

        log('Resolving %s.depend' % self.name, 7)
        if config.preview is False:
            for j in self.depend:
                j.results.wait()
        log('%s.depend resolved' % self.name, 7)

        jobmanager._register(self)

        log('Starting %s.prerun()' % self.name, 5)
        self.prerun()
        log('%s.prerun() finished' % self.name, 5)

        for i in reversed(self.default_settings):
            self.settings.soft_update(i)

        prev = jobmanager._check_hash(self)
        if prev is not None:
            try:
                prev.results._copy_to(self.results)
                self.status = 'copied'
            except ResultsError as re:
                log('Copying results of %s failed because of the following error: %s' % (prev.name, str(re)), 1)
                self.status = prev.status
            if self.settings.pickle:
                self.pickle()
            self.results.finished.set()
            self.results.done.set()
            if self.parent:
                self.parent._notify()
        else:
            self.status = 'running'
            log('Starting %s._get_ready()' % self.name, 7)
            self._get_ready()
            log('%s._get_ready() finished' % self.name, 7)

        log('%s._prepare() finished' % self.name, 7)
        return prev is None



    def _get_ready(self):
        """Get ready for :meth:`~Job._execute`. This is the last step before :meth:`~Job._execute` is called. Abstract method."""
        raise PlamsError('Trying to run an abstract method Job._get_ready()')

    def _execute(self, jobrunner):
        """Execute the job. Abstract method."""
        raise PlamsError('Trying to run an abstract method Job._execute()')



    def _finalize(self):
        """Gather the results of job execution and organize them. This method collects steps 9-12 from :ref:`job-life-cycle`. Should not be overridden."""
        log('Starting %s._finalize()' % self.name, 7)

        if config.preview is False:
            log('Collecting results of %s' % self.name, 7)
            self.results.collect()
            self.results.finished.set()
            if self.status != 'crashed':
                self.status = 'finished'
                if self.check():
                    log('%s.check() success. Cleaning results with keep = %s' % (self.name, self.settings.keep), 7)
                    self.results._clean(self.settings.keep)
                    log('Starting %s.postrun()' % self.name, 5)
                    self.postrun()
                    log('%s.postrun() finished' % self.name, 5)
                    self.status = 'successful'
                    log('Pickling %s' % self.name, 7)
                    if self.settings.pickle:
                        self.pickle()
                else:
                    log('%s.check() failed' % self.name, 7)
                    self.status = 'failed'
            self.results.done.set()
        else:
            self.status = 'preview'

        if self.parent:
            self.parent._notify()

        log('%s._finalize() finished' % self.name, 7)
        log("Job %s finished with status '%s' "% (self.name, self.status), 1)


#===================================================================================================
#===================================================================================================
#===================================================================================================


class SingleJob(Job):
    """Abstract class representing a job consisting of a single execution of some external binary (or arbitrary shell script in general).

    In addition to constructor arguments and attributes defined by |Job|, the constructor of this class accepts the keyword argument ``molecule`` that should be a |Molecule| instance.

    Class attribute ``_filenames`` defines default names for input, output, runscript and error files. If you wish to override this attribute it should be a dictionary with string keys ``'inp'``, ``'out'``, ``'run'``, ``'err'``. The value for each key should be a string describing corresponding file's name. Shortcut ``$JN`` can be used for job's name. The default value is defined in the following way::

        >>> _filenames = {'inp':'$JN.in', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}

    This class defines no new methods that could be directly called in your script. Methods that can and should be overridden are :meth:`~SingleJob.get_input` and :meth:`~SingleJob.get_runscript`.

    """
    _filenames = {'inp':'$JN.in', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}

    def __init__(self, molecule=None, **kwargs):
        Job.__init__(self, **kwargs)
        self.molecule = molecule



    def _filename(self, t):
        """Return filename for file of type *t*. *t* can be any key from ``_filenames`` dictionary. ``$JN`` is replaced with job name in returned string."""
        return self.__class__._filenames[t].replace('$JN', self.name)



    def get_input(self):
        """Generate the input file. Abstract method.

        This method should return a single string with full content of the input file. It should process information stored in ``input`` branch of job's settings and in ``molecule`` attribute.
        """
        raise PlamsError('Trying to run an abstract method SingleJob._get_input()')

    def get_runscript(self):
        """Generate runscript. Abstract method.

        This method should return a single string with runscript contents. It can process information stored in ``runscript`` branch of job's settings. In general the full runscript has the following form::

            [first line defined by job.settings.runscript.shebang]

            [contents of job.settings.runscript.pre, if any]

            [value returned by get_runscript()]

            [contents of job.settings.runscript.post, if any]

        When overridden, this method should pay attention to ``.runscript.stdout_redirect`` key in job's ``settings``.
        """
        raise PlamsError('Trying to run an abstract method SingleJob._get_runscript()')



    def hash(self):
        """Calculate unique hash of this instance.

        The behavior of this method is adjusted by the value of ``hashing`` key in |JobManager| settings. If no |JobManager| is yet associated with this job, default setting from ``config.jobmanager.hashing`` is used.

        Currently supported values for ``hashing`` are:
            *   ``False`` or ``None`` -- returns ``None`` and disables |RPM|.
            *   ``input`` -- returns SHA256 hash of the input file.
            *   ``runscript`` -- returns SHA256 hash of the runscript.
            *   ``input+runscript`` -- returns SHA256 hash of the concatenation of input and runscript.
        """

        if self.jobmanager:
            mode = self.jobmanager.settings.hashing
        else:
            mode = config.jobmanager.hashing

        if not mode:
            return None
        h = hashlib.sha256()
        if mode == 'input':
            h.update(self.get_input().encode())
        elif mode == 'runscript':
            h.update(self._full_runscript().encode())
        elif mode == 'input+runscript':
            h.update(self.get_input().encode())
            h.update(self._full_runscript().encode())
        else:
            raise PlamsError('Unsupported hashing method: ' + str(mode))
        return h.hexdigest()


    def _full_runscript(self):
        """Generate full runscript, including shebang line and contents of ``pre`` and ``post``, if any.

        .. technical::

            In practice this method is just a wrapper around :meth:`~SingleJob.get_runscript`.
        """
        ret = self.settings.runscript.shebang +'\n\n'
        if 'pre' in self.settings.runscript:
            ret += self.settings.runscript.pre+'\n\n'
        ret += self.get_runscript()
        if 'post' in self.settings.runscript:
            ret += self.settings.runscript.post+'\n\n'
        return ret


    def _get_ready(self):
        """Generate input and runscript files in the job folder. Methods :meth:`get_input` and :meth:`get_runscript` are used for that purpose."""
        inpfile = opj(self.path, self._filename('inp'))
        runfile = opj(self.path, self._filename('run'))

        with open(inpfile, 'w') as inp:
            inp.write(self.get_input())

        with open(runfile, 'w') as run:
            run.write(self._full_runscript())

        os.chmod(runfile, os.stat(runfile).st_mode | stat.S_IEXEC)


    def _execute(self, jobrunner):
        """Execute previously created runscript using *jobrunner*.

        The method :meth:`~scm.plams.jobrunner.JobRunner.call` of *jobrunner* is used. Working directory is ``self.path``. ``self.settings.run`` is passed as ``runflags`` argument.

        If preview mode is on, this method does nothing.
        """
        log('Starting %s._execute()' % self.name, 7)
        if config.preview is False:
            o = self._filename('out') if not self.settings.runscript.stdout_redirect else None
            retcode = jobrunner.call(runscript=self._filename('run'), workdir=self.path, out=o, err=self._filename('err'), runflags=self.settings.run)
            if retcode != 0:
                log('WARNING: Job %s finished with nonzero return code' % self.name, 1)
                self.status = 'crashed'
        log('%s._execute() finished' % self.name, 7)


#===================================================================================================
#===================================================================================================
#===================================================================================================


class MultiJob(Job):
    """Concrete class representing a job that is a container for other jobs.

    In addition to constructor arguments and attributes defined by |Job|, the constructor of this class accepts two keyword arguments:
        *   ``children`` -- should be a list (or other iterable container) containing children jobs.
        *   ``childrunner`` -- by default all the children jobs are run using the same |JobRunner| as the parent job. If you wish to use a different |JobRunner| for children, you can pass it using this argument.

    Values passed as ``children`` and ``childrunner`` are stored as instance attributes and can be adjusted later, but before the |run| method is called.

    This class defines no new methods that could be directly called in your script.

    When executed, a multijob runs all its children using the same |run| arguments. If you need to specify different run flags for children you can do it by manually setting them in children job |Settings|::

        >>> childjob.settings.run.arg = 'value'

    Since ``run`` branch of settings gets soft-updated by run flags, value set this way is not overwritten by parent job.

    Job folder of a multijob gets cleaned independently of its children. See |cleaning| for details.
    """
    def __init__(self, children=None, childrunner=None, **kwargs):
        Job.__init__(self, **kwargs)
        self.children = [] if children is None else children
        self.childrunner = childrunner
        self._active_children = 0
        self._lock =  threading.Lock()


    def new_children(self):
        """Generate new children jobs.

        This method is useful when some of children jobs are not known beforehand and need to be generated based on other children jobs, like for example in any kind of self-consistent procedure.

        The goal of this method is to produce new portion of children jobs. Newly created jobs **have to** be manually added to ``self.children`` and, besides that, returned as a list by this method. No adjustment of newly created jobs' ``parent`` attribute is needed. This method **cannot** modify ``_active_children`` attribute.

        The method defined here is a default template, returning an empty list, which means no new children jobs are generated and the entire execution of the parent job consists only of running jobs initially found in ``self.children``. To modify this behavior you can override this method in |MultiJob| subclass or use one of |binding_decorators|, just like with :ref:`prerun-postrun`.
        """

        return []


    def hash(self):
        """Hashing for multijobs is disabled by default. Return ``None``."""
        return None


    def check(self):
        """Check if the calculation was successful. Returns ``True`` if every children job has its ``status`` attribute set to ``'successful'``.
        """
        return all([child.status in ['successful', 'copied'] for child in self])


    def _get_ready(self):
        """Get ready for :meth:`~MultiJob._execute`. Count children jobs and set their ``parent`` attribute."""
        self._active_children = len(self.children)
        for child in self:
            child.parent = self


    def __iter__(self):
        """Iterate through ``children``. If it is a dictionary, iterate through its values."""
        if isinstance(self.children, dict):
            return iter(self.children.values())
        return iter(self.children)


    def _notify(self):
        """Notify this job that one of its children has finished.

        Decrement ``_active_children`` by one. Use ``_lock`` to ensure thread safety.
        """
        with self._lock:
            self._active_children -= 1


    def _execute(self, jobrunner):
        """Run all children from ``children``. Then use :meth:`~MultiJob.new_children` and run all jobs produced by it. Repeat this procedure until :meth:`~MultiJob.new_children` returns an empty list. Wait for all started jobs to finish."""
        log('Starting %s._execute()' % self.name, 7)
        jr = self.childrunner or jobrunner
        for child in self:
            child.run(jobrunner=jr, jobmanager=self.jobmanager, **self.settings.run)
        new = self.new_children()
        while new:
            with self._lock:
                self._active_children += len(new)
            for child in new:
                child.parent = self
                child.run(jobrunner=jr, jobmanager=self.jobmanager, **self.settings.run)
            new = self.new_children()
        while self._active_children > 0:
            time.sleep(config.sleepstep)
        log('%s._execute() finished' % self.name, 7)
