import os
import stat
import threading
import time

try:
    import dill as pickle
except ImportError:
    import pickle

from os.path import join as opj

from .errors import FileError, JobError, PlamsError, ResultsError
from .functions import config, log
from .private import sha256
from .results import Results
from .settings import Settings
from ..mol.molecule import Molecule

__all__ = ['SingleJob', 'MultiJob']


class Job:
    """General abstract class for all kind of computational tasks.

    Methods common for all kinds of jobs are gathered here. Instances of |Job| should never be created. It should not be subclassed either. If you wish to define a new type of job please subclass either |SingleJob| or |MultiJob|.

    Methods that are meant to be explicitly called by the user are |run| and occasionally :meth:`~Job.pickle`. In most cases |pickling| is done automatically, but if for some reason you wish to do it manually, you can use :meth:`~Job.pickle` method.

    Methods that can be safely overridden in subclasses are:

    *   :meth:`~Job.check`
    *   :meth:`~Job.hash` (see |RPM|)
    *   |prerun| and |postrun| (see :ref:`prerun-postrun`)

    All other methods should remain unchanged.

    Class attribute ``_result_type`` defines the type of results associated with this job. It should point to a class and it **must** be a |Results| subclass.

    Every job instance has the following attributes:

    Attributes adjusted automatically that should not be changed by the user:

    *   ``status`` -- current status of the job. Possible values: *created*, *started*, *registered*, *running*, *finished*, *crashed*, *failed*, *successful*, *copied*, *preview*.
    *   ``results`` -- a reference to a results instance. An empty instance of the type stored in ``_result_type`` is created when the job object is created.
    *   ``path`` -- an absolute path to the job folder.
    *   ``jobmanager`` -- a job manager associated with this job.
    *   ``parent`` -- a pointer to the parent job if this job is a child job of some |MultiJob|. ``None`` otherwise.

    Attributes that can be modified, but only before |run| is called:

    *   ``name`` -- the name of the job.
    *   ``settings`` -- settings of the job.
    *   ``default_settings`` -- see :ref:`default-settings`.
    *   ``depend`` -- a list of explicit dependencies.
    *   ``_dont_pickle`` -- additional list of this instance's attributes that will be removed before pickling. See |pickling| for details.

    """

    _result_type = Results

    def __init__(self, name='plamsjob', settings=None, depend=None):
        if os.path.sep in name:
            raise PlamsError('Job name cannot contain {}'.format(os.path.sep))
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


    #=======================================================================


    def run(self, jobrunner=None, jobmanager=None, **kwargs):
        """Run the job using *jobmanager* and *jobrunner* (or defaults, if ``None``). Other keyword arguments (*\*\*kwargs*) are stored in ``run`` branch of job's settings. Returned value is the |Results| instance associated with this job.

        .. warning::

            This method should **not** be overridden.

        .. technical::

            This method does not do too much by itself. After some initial preparation it passes control to the job runner, which decides if a new thread should be started for this job. The role of the job runner is to execute three methods that make the full job life cycle: :meth:`~Job._prepare`, :meth:`~Job._execute` and :meth:`~Job._finalize`. During :meth:`~Job._execute` the job runner is called once again to execute the runscript (only in case of |SingleJob|).
        """
        if self.status != 'created':
            raise JobError('Trying to run previously started job {}'.format(self.name))

        self.status = 'started'
        self._log_status(1)

        self.settings.run.soft_update(Settings(kwargs))

        if jobrunner is None:
            if 'default_jobrunner' in config:
                jobrunner = config.default_jobrunner
            else:
                raise PlamsError('No default jobrunner found. This probably means that PLAMS init() was not called.')
        if jobmanager is None:
            if 'default_jobmanager' in config:
                jobmanager = config.default_jobmanager
            else:
                raise PlamsError('No default jobmanager found. This probably means that PLAMS init() was not called.')

        jobrunner._run_job(self, jobmanager)
        return self.results


    def pickle(self, filename=None):
        """Pickle this instance and save to a file indicated by *filename*. If ``None``, save to ``[jobname].dill`` in the job folder."""
        filename = filename or opj(self.path, self.name+'.dill')
        with open(filename, 'wb') as f:
            try:
                pickle.dump(self, f, -1)
            except:
                log('Pickling of {} failed'.format(self.name), 1)


    def ok(self, strict=True):
        """Check if the execution of this instance was successful. If needed, wait for the job to finish and then check if the status is *successful* (or *copied*).

            If this method is called before job's |run| method, a warning is logged and the returned value is ``False``. The most likely cause of such a behavior is simply forgetting about |run| call. However, in complicated workflows executed in parallel, it can sometimes naturally happen that one thread is ahead of others and calls :meth:`~Job.ok` before some other thread has a chance to call |run|. If you're experiencing that kind of problems, please consider using ``strict=False`` to skip the |run| check.  But keep in mind that skipping that check will deadlock the current thread if |run| never gets called.

        """
        if strict and self.status == 'created': #first thing run() does is changing the status to 'started'
            log('Job {} WARNING: ok() method was called before run(). Returned value is False. Please check the documentation'.format(self.name), 3)
            return False
        self.results.wait()
        return self.status in ['successful', 'copied']


    def check(self):
        """Check if the execution of this instance was successful. Abstract method meant for internal use."""
        raise PlamsError('Trying to run an abstract method Job.check()')


    def hash(self):
        """Calculate the hash of this instance. Abstract method meant for internal use."""
        raise PlamsError('Trying to run an abstract method Job.hash()')


    def prerun(self):
        """Actions to take before the actual job execution.

        This method is initially empty, it can be defined in subclasses or directly added to either the whole class or a single instance using |binding_decorators|.
        """
        pass


    def postrun(self):
        """Actions to take just after the actual job execution.

        This method is initially empty, it can be defined in subclasses or directly added to either the whole class or a single instance using |binding_decorators|.
        """
        pass


    #=======================================================================


    def _prepare(self, jobmanager):
        """Prepare the job for execution. This method collects steps 1-7 from :ref:`job-life-cycle`. Should not be overridden. Returned value indicates if job execution should continue (|RPM| did not find this job as previously run)."""

        log('Starting {}._prepare()'.format(self.name), 7)

        log('Resolving {}.depend'.format(self.name), 7)
        if config.preview is False:
            for j in self.depend:
                j.results.wait()
        log('{}.depend resolved'.format(self.name), 7)

        jobmanager._register(self)

        log('Starting {}.prerun()'.format(self.name), 5)
        self.prerun()
        log('{}.prerun() finished'.format(self.name), 5)

        for i in reversed(self.default_settings):
            self.settings.soft_update(i)

        prev = jobmanager._check_hash(self)
        if prev is not None:
            try:
                prev.results._copy_to(self.results)
                self.status = 'copied'
            except ResultsError as re:
                log('Copying results of {} failed because of the following error: {}'.format(prev.name, str(re)), 1)
                self.status = prev.status
            if self.settings.pickle:
                self.pickle()
            self.results.finished.set()
            self.results.done.set()
            if self.parent and self in self.parent:
                self.parent._notify()
        else:
            self.status = 'running'
            log('Starting {}._get_ready()'.format(self.name), 7)
            self._get_ready()
            log('{}._get_ready() finished'.format(self.name), 7)

        log('{}._prepare() finished'.format(self.name), 7)
        self._log_status(3)
        return prev is None


    def _get_ready(self):
        """Get ready for :meth:`~Job._execute`. This is the last step before :meth:`~Job._execute` is called. Abstract method."""
        raise PlamsError('Trying to run an abstract method Job._get_ready()')


    def _execute(self, jobrunner):
        """Execute the job. Abstract method."""
        raise PlamsError('Trying to run an abstract method Job._execute()')


    def _finalize(self):
        """Gather the results of the job execution and organize them. This method collects steps 9-12 from :ref:`job-life-cycle`. Should not be overridden."""
        log('Starting {}._finalize()'.format(self.name), 7)

        if config.preview is False:
            log('Collecting results of {}'.format(self.name), 7)
            self.results.collect()
            self.results.finished.set()
            if self.status != 'crashed':
                self.status = 'finished'
                self._log_status(3)
                if self.check():
                    log('{}.check() success. Cleaning results with keep = {}'.format(self.name, self.settings.keep), 7)
                    self.results._clean(self.settings.keep)
                    log('Starting {}.postrun()'.format(self.name), 5)
                    self.postrun()
                    log('{}.postrun() finished'.format(self.name), 5)
                    self.status = 'successful'
                    log('Pickling {}'.format(self.name), 7)
                    if self.settings.pickle:
                        self.pickle()
                else:
                    log('{}.check() failed'.format(self.name), 7)
                    self.status = 'failed'
        else:
            self.status = 'preview'
            self.results.finished.set()
        self.results.done.set()

        if self.parent and self in self.parent:
            self.parent._notify()

        log('{}._finalize() finished'.format(self.name), 7)
        self._log_status(1)


    def __getstate__(self):
        """Prepare this job instance for pickling.

        Attributes ``jobmanager``, ``parent``, ``default_settings`` and ``_lock`` are removed, as well as all attributes listed in ``self._dont_pickle``.
        """
        remove = ['jobmanager', 'parent', 'default_settings', '_lock'] + self._dont_pickle
        return {k:v for k,v in self.__dict__.items() if k not in remove}


    def _log_status(self, level):
        """Log the status of this instance on a chosen log *level*. The message is uppercased to clearly stand out among other log entries."""
        log('JOB {} {}'.format(self._full_name(), self.status.upper()), level)


    def _full_name(self):
        if self.parent:
            return '/'.join([self.parent._full_name(), self.name])
        return self.name


#===========================================================================
#===========================================================================
#===========================================================================


class SingleJob(Job):
    """Abstract class representing a job consisting of a single execution of some external binary (or arbitrary shell script in general).

    In addition to constructor arguments and attributes defined by |Job|, the constructor of this class accepts the keyword argument ``molecule`` that should be a |Molecule| instance. The constructor creates a copy of the supplied |Molecule| and stores it as the ``molecule`` attribute.

    Class attribute ``_filenames`` defines default names for input, output, runscript and error files. If you wish to override this attribute it should be a dictionary with string keys ``'inp'``, ``'out'``, ``'run'``, ``'err'``. The value for each key should be a string describing corresponding file's name. Shortcut ``$JN`` can be used for job's name. The default value is defined in the following way::

        _filenames = {'inp':'$JN.in', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}

    This class defines no new methods that could be directly called in your script. Methods that can and should be overridden are |get_input| and |get_runscript|.

    """
    _filenames = {'inp':'$JN.in', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}

    def __init__(self, molecule=None, **kwargs):
        Job.__init__(self, **kwargs)
        self.molecule = molecule.copy() if isinstance(molecule, Molecule) else molecule


    def get_input(self):
        """Generate the input file. Abstract method.

        This method should return a single string with the full content of the input file. It should process information stored in the ``input`` branch of job settings and in the ``molecule`` attribute.
        """
        raise PlamsError('Trying to run an abstract method SingleJob._get_input()')


    def get_runscript(self):
        """Generate the runscript. Abstract method.

        This method should return a single string with the runscript contents. It can process information stored in ``runscript`` branch of job  settings. In general the full runscript has the following form::

            [first line defined by job.settings.runscript.shebang]

            [contents of job.settings.runscript.pre, when present]

            [value returned by get_runscript()]

            [contents of job.settings.runscript.post, when present]

        When overridden, this method should pay attention to ``.runscript.stdout_redirect`` key in job's ``settings``.
        """
        raise PlamsError('Trying to run an abstract method SingleJob._get_runscript()')


    def hash_input(self):
        """Calculate SHA256 hash of the input file."""
        return sha256(self.get_input())


    def hash_runscript(self):
        """Calculate SHA256 hash of the runscript."""
        return sha256(self.full_runscript())


    def hash(self):
        """Calculate unique hash of this instance.

        The behavior of this method is adjusted by the value of ``hashing`` key in |JobManager| settings. If no |JobManager| is yet associated with this job, default setting from ``config.jobmanager.hashing`` is used.

        Methods :meth:`~SingleJob.hash_input` and :meth:`~SingleJob.hash_runscript` are used to obtain hashes of, respectively, input and runscript.

        Currently supported values for ``hashing`` are:

        *   ``False`` or ``None`` -- returns ``None`` and disables |RPM|.
        *   ``input`` -- returns the hash of the input file.
        *   ``runscript`` -- returns the hash of the runscript.
        *   ``input+runscript`` -- returns SHA256 hash of the concatenation of **hashes** of input and runscript.
        """
        if self.jobmanager:
            mode = self.jobmanager.settings.hashing
        else:
            mode = config.jobmanager.hashing

        if not mode:
            return None
        if mode == 'input':
            return self.hash_input()
        elif mode == 'runscript':
            return self.hash_runscript()
        elif mode == 'input+runscript':
            return sha256(self.hash_input() + self.hash_runscript())
        else:
            raise PlamsError('Unsupported hashing method: {}'.format(mode))


    def check(self):
        """Check if the calculation was successful.

        This method can be overridden in concrete subclasses of |SingleJob|. It should return a boolean value. The definition here serves as a default, to prevent crashing if a subclass does not define its own :meth:`~scm.plams.core.basejob.SingleJob.check`. It always returns ``True``.

        .. warning::

            This method is meant for internal usage and **should not** be explicitly called in your script (but it can be overridden in subclasses). Manually calling :meth:`~scm.plams.core.basejob.SingleJob.check` is not thread safe. For a thread safe function to evaluate the state of your job please use :meth:`~scm.plams.core.basejob.Job.ok`

        """
        return True


    def full_runscript(self):
        """Generate the full runscript, including shebang line and contents of ``pre`` and ``post``, if any. In practice this method is just a simple wrapper around |get_runscript|.
        """
        ret = self.settings.runscript.shebang +'\n\n'
        if 'pre' in self.settings.runscript:
            ret += self.settings.runscript.pre+'\n\n'
        ret += self.get_runscript()
        if 'post' in self.settings.runscript:
            ret += self.settings.runscript.post+'\n\n'
        return ret


    def _get_ready(self):
        """Create input and runscript files in the job folder. Methods |get_input| and :meth:`full_runscript` are used for that purpose. Filenames correspond to entries in the `_filenames` attribute"""
        inpfile = opj(self.path, self._filename('inp'))
        runfile = opj(self.path, self._filename('run'))

        with open(inpfile, 'w') as inp:
            inp.write(self.get_input())

        with open(runfile, 'w') as run:
            run.write(self.full_runscript())

        os.chmod(runfile, os.stat(runfile).st_mode | stat.S_IEXEC)


    def _execute(self, jobrunner):
        """Execute previously created runscript using *jobrunner*.

        The method :meth:`~scm.plams.core.jobrunner.JobRunner.call` of *jobrunner* is used. Working directory is ``self.path``. ``self.settings.run`` is passed as ``runflags`` argument.

        If preview mode is on, this method does nothing.
        """
        log('Starting {}._execute()'.format(self.name), 7)
        if config.preview is False:
            o = self._filename('out') if not self.settings.runscript.stdout_redirect else None
            retcode = jobrunner.call(runscript=self._filename('run'), workdir=self.path, out=o, err=self._filename('err'), runflags=self.settings.run)
            if retcode != 0:
                log('WARNING: Job {} finished with nonzero return code'.format(self.name), 3)
                self.status = 'crashed'
        log('{}._execute() finished'.format(self.name), 7)


    def _filename(self, t):
        """Return filename for file of type *t*. *t* can be any key from ``_filenames`` dictionary. ``$JN`` is replaced with job name in the returned string."""
        return self._filenames[t].replace('$JN', self.name)


    @classmethod
    def load_external(cls, path, settings=None, molecule=None, finalize=False):
        """Load an external job from *path*.

        In this context an "external job" is an execution of some external binary that was not managed by PLAMS, and hence does not have a ``.dill`` file. It can also be used in situations where the execution was started with PLAMS, but the Python process was terminated before the execution finished, resulting in steps 9-12 of :ref:`job-life-cycle` not happening.

        All the files produced by your computation should be placed in one folder and *path* should be the path to this folder. The name of the folder is used as a job name. Input, output, error and runscript files, if present, should have names defined in ``_filenames`` class attribute (usually ``[jobname].in``, ``[jobname].out``, ``[jobname].err`` and ``[jobname].run``). It is not required to supply all these files, but in most cases one would like to use at least the output file, in order to use methods like :meth:`~scm.plams.core.results.Results.grep_output` or :meth:`~scm.plams.core.results.Results.get_output_chunk`.

        This method is a class method, so it is called via class object and it returns an instance of that class::

            >>> j = SingleJob.load_external(path='some/path/jobname')
            >>> type(j)
            scm.plams.core.basejob.SingleJob
            >>> a = AMSJob.load_external(path='some/path/jobname')
            >>> type(a)
            scm.plams.interfaces.adfsuite.ams.AMSJob

        You can supply |Settings| and |Molecule| instances as *settings* and *molecule* parameters, they will end up attached to the returned job instance. If you don't do this, PLAMS will try to recreate them automatically using methods :meth:`~scm.plams.core.results.Results.recreate_settings` and :meth:`~scm.plams.core.results.Results.recreate_molecule` of the corresponding |Results| subclass. If no |Settings| instance is obtained in either way, the defaults from ``config.job`` are copied.

        You can set the *finalize* parameter to ``True`` if you wish to run the whole :meth:`~Job._finalize` on the newly created job. In that case PLAMS will perform the usual :meth:`~Job.check` to determine the job status (*successful* or *failed*), followed by cleaning of the job folder (|cleaning|), |postrun| and pickling (|pickling|). If *finalize* is ``False``, the status of the returned job is *copied*.
        """

        if not os.path.isdir(path):
            raise FileError('Path {} does not exist, cannot load from it.'.format(path))

        path = os.path.abspath(path)
        jobname = os.path.basename(path)

        job = cls(name=jobname)
        job.path = path
        job.status = 'copied'
        job.results.collect()

        job._filenames = {}
        for t,name in cls._filenames.items():
            fullname = name.replace('$JN', job.name)
            if fullname in job.results.files:
                job._filenames[t] = name
            else:
                log('The default {} file {} not present in {}'.format(t, fullname, job.path), 5)

        job.settings  = settings or job.results.recreate_settings() or config.job.copy()
        job.molecule  = molecule or job.results.recreate_molecule()

        if finalize:
            job._finalize()

        return job




#===========================================================================
#===========================================================================
#===========================================================================


class MultiJob(Job):
    """Concrete class representing a job that is a container for other jobs.

    In addition to constructor arguments and attributes defined by |Job|, the constructor of this class accepts two keyword arguments:

    *   ``children`` -- iterable container with children jobs (usually a list or a dictionary).
    *   ``childrunner`` -- by default all the children jobs are run using the same |JobRunner| as the parent job. If you wish to use a different |JobRunner| for children, you can pass it using ``childrunner``.

    Values passed as ``children`` and ``childrunner`` are stored as instance attributes and can be adjusted later, but before the |run| method is called.

    This class defines no new methods that could be directly called in your script.

    When executed, a multijob runs all its children using the same |JobManager| and |JobRunner| (unless ``childrunner`` is set). If you need to specify different run flags for children you can do it by manually setting them in child job settings::

        childjob = myjob.children[0]
        childjob.settings.run.arg = 'value'

    Since the ``run`` branch of settings gets soft-updated with run flags, the value set this way is not overwritten by parent job.

    The job folder of a multijob gets cleaned independently of its children. See |cleaning| for details.

    Private attributes ``_active_children`` and ``_lock`` are essential for proper parallel execution. Please do not modify them.
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

        The goal of this method is to produce a new portion of children jobs. Newly created jobs should be returned in a container compatible with ``self.children`` (e.g. list for list, dict for dict). No adjustment of newly created jobs' ``parent`` attribute is needed. This method **cannot** modify ``_active_children`` attribute.

        The method defined here is a default template, returning ``None``, which means no new children jobs are generated and the entire execution of the parent job consists only of running jobs initially found in ``self.children``. To modify this behavior you can override this method in a |MultiJob| subclass or you can use one of |binding_decorators|, just like with :ref:`prerun-postrun`.
        """
        return None


    def hash(self):
        """Hashing for multijobs is disabled by default. Returns ``None``."""
        return None


    def check(self):
        """Check if the execution of this instance was successful, by calling :meth:`Job.ok` of all the children jobs.
        """
        return all([child.ok() for child in self])


    def other_jobs(self):
        """Iterate through other jobs that belong to this |MultiJob|, but are not in ``children``.

        Sometimes |prerun| or |postrun| methods create and run some small jobs that don't end up in ``children`` collection, but are still considered a part of a |MultiJob| instance (their ``parent`` atribute points to the |MultiJob| and their working folder is inside MultiJob's working folder). This method provides an iterator that goes through all such jobs.

        Each attribute of this |MultiJob| that is of type |Job| and has it's parent pointing to this |MultiJob| is returned, in a random order.
        """
        for attr in self.__dict__.values():
            if isinstance(attr, Job) and ((hasattr(attr, 'parent') and attr.parent == self) or not hasattr(attr, 'parent')):
                yield attr


    def remove_child(self, job):
        """Remove *job* from children."""

        rm = None
        for i, j in (self.children.items() if isinstance(self.children, dict) else enumerate(self.children)):
            if j == job:
                rm = i
                break
        if rm is not None:
            del self.children[rm]


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
        log('Starting {}._execute()'.format(self.name), 7)
        jr = self.childrunner or jobrunner

        for child in self:
            child.run(jobrunner=jr, jobmanager=self.jobmanager, **self.settings.run)

        new = self.new_children()
        while new:
            with self._lock:
                self._active_children += len(new)

            if isinstance(new, dict) and isinstance(self.children, dict):
                self.children.update(new)
                it = new.values()
            elif isinstance(new, list) and isinstance(self.children, list):
                self.children += new
                it = new
            else:
                raise JobError("ERROR in job {}: 'new_children' returned a value incompatible with 'children'".format(self.name))

            for child in it:
                child.parent = self
                child.run(jobrunner=jr, jobmanager=self.jobmanager, **self.settings.run)

            new = self.new_children()

        while self._active_children > 0:
            time.sleep(config.sleepstep)
        log('{}._execute() finished'.format(self.name), 7)
