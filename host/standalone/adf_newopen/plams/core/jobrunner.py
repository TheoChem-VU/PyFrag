from __future__ import unicode_literals
from six import add_metaclass

import os
import functools
import threading
import time
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

from os.path import join as opj

from .common import log, string
from .errors import PlamsError
from .settings import Settings
from .basejob import SingleJob


__all__ = ['JobRunner', 'GridRunner']


def _in_thread(func):
    """Decorator for an instance method. If ``parallel`` attribute of given instance is ``True``, run decorated method in a separate :class:`~threading.Thread`. This thread is usually a daemon thread, the decision is based on ``config.daemon_threads`` entry."""
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.parallel:
            t = threading.Thread(name='plamsthread', target=func, args=(self,)+args, kwargs=kwargs)
            t.daemon = config.daemon_threads
            t.start()
        else:
            func(self, *args, **kwargs)
    return wrapper

def _limit(func):
    """Decorator for an instance method. If ``semaphore`` attribute of given instance is not ``None``, use this attribute to wrap decorated method via :ref:`with<with-locks>` statement."""
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.semaphore:
            with self.semaphore:
                return func(self, *args, **kwargs)
        else:
            return func(self, *args, **kwargs)
    return wrapper

class _MetaRunner(type):
    """Metaclass for |JobRunner|. Wraps :meth:`~scm.plams.jobrunner.JobRunner.call` with :func:`_limit` decorator."""
    def __new__(meta, name, bases, dct):
        dct['call'] = _limit(dct['call'])
        return type.__new__(meta, name, bases, dct)

@add_metaclass(_MetaRunner)
class JobRunner(object):
    """Class representing local job runner.

    The goal of |JobRunner| is to take care of two important things -- parallelization and runscript execution:
        *   When the method |run| of any |Job| instance is executed, this method, after some preparations, passes control to a |JobRunner| instance. This |JobRunner| instance decides if a separate thread should be spawned for this job or if the execution should proceed in the main thread. This decision is based on ``parallel`` attribute which can be set on |JobRunner| creation. There are no separate classes for serial and parallel job runner, both cases are covered by |JobRunner| depending on one bool parameter.
        *   If the executed job is an instance of |SingleJob|, it creates a shell script (called runscript) which contains most of the actual computational work (usually it is just an execution of some external binary). The runscript is then submitted to a |JobRunner| instance using its method :meth:`call`. This method executes the runscript in a separate subprocess and takes care of setting proper working directory, output and error stream handling etc.

    The number of simultaneously running :meth:`call` methods can be limited using *maxjobs* parameter. If *maxjobs* is 0, no limit is enforced. If *parallel* is ``False``, *maxjobs* is ignored. If *parallel* is ``True`` and *maxjobs* is a positive integer, a :func:`BoundedSemaphore<threading.BoundedSemaphore>` of that size is used to limit the number of concurrently running :meth:`call` methods.

    A |JobRunner| instance can be passed to |run| with a keyword argument ``jobrunner``. If this argument is omitted, the instance stored in ``config.default_jobrunner`` is used.
    """

    def __init__ (self, parallel=False, maxjobs=0):
        self.parallel = parallel
        self.semaphore = threading.BoundedSemaphore(maxjobs) if maxjobs else None


    def call(self, runscript, workdir, out, err, **kwargs):
        """call(runscript, workdir, out, err, **kwargs)
        Execute the *runscript* in the folder *workdir*. Redirect output and error streams to *out* and *err*, respectively.

        Arguments mentioned above should be strings containing paths to corresponding files or folders

        Other keyword arguments are ignored here but they can be useful in |JobRunner| subclasses (see :meth:`GridRunner.call`).

        Returns integer value indicating the exit code returned by execution of *runscript*.

        This method can be safely overridden in |JobRunner| subclasses. For example, in |GridRunner| it submits the runscript to a job scheduler instead of executing it locally.

        .. note::
            This method is used automatically during |run| and should never be explicitly called in your script.
        """
        log('Executing %s' % runscript, 5)
        with open(opj(workdir, err), 'w') as e:
            kwargs = {'cwd':workdir, 'stderr':e}
            if os.name == 'posix':
                command = './'+runscript
                kwargs['close_fds'] = True
            else:
                command = ['sh', runscript]
            if out is not None:
                kwargs['stdout'] = open(opj(workdir, out), 'w')
            retcode = subprocess.call(command, **kwargs)
            if out is not None:
                kwargs['stdout'].close()
        log('Execution of %s finished with returncode %i' % (runscript, retcode), 5)
        return retcode

    @_in_thread
    def _run_job(self, job, jobmanager):
        """_run_job(job, jobmanager)
        This method aggregates these parts of |run| that are supposed to be run in a separate thread in case of parallel job execution. It is wrapped with :func:`_in_thread` decorator.

        This method should not be overridden.
        """
        if job._prepare(jobmanager):
            job._execute(self)
            job._finalize()

#===================================================================================================
#===================================================================================================
#===================================================================================================


class GridRunner(JobRunner):
    """Subclass of |JobRunner| that submits the runscript to a job scheduler instead of executing it locally. Besides two new keyword arguments (*grid* and *sleepstep*) and different :meth:`call` method it behaves and is meant to be used just like regular |JobRunner|.

    There are many different job schedulers that are popular and widely used nowadays (for example TORQUE, SLURM, OGE). Usually they use different commands for submitting jobs or checking queue status. This class tries to build a common and flexible interface for all those tools. The idea is that commands used to communicate with job scheduler are not rigidly hard-coded but dynamically taken from a |Settings| instance instead. Thanks to that user has almost full control over the behavior of |GridRunner|.

    So the behavior of |GridRunner| is determined by the contents of |Settings| instance stored in its ``settings`` attribute. This |Settings| instance can be manually supplied by the user or taken from a collection of predefined behaviors stored as branches of ``config.gridrunner``. The adjustment is done via *grid* parameter that should be either string or |Settings|. If string, it has to be a key occurring in ``config.gridrunner`` (or ``'auto'`` for autodetection). For example, if ``grid='slurm'`` is passed, ``config.gridrunner.slurm`` is linked as ``settings``. If *grid* is ``'auto'``, entries in ``config.gridrunner`` are tested one by one and the first one that works (its submit command is present on your system) is chosen. When a |Settings| instance is passed it gets plugged directly as ``settings``.

    Currently two predefined job schedulers are available (see ``plams_defaults.py``): ``slurm`` for SLURM and ``pbs`` for job schedulers following PBS syntax (PBS, TORQUE, Oracle Grid Engine etc.).

    The |Settings| instance used for |GridRunner| should have the following structure:
        *   ``.output`` -- flag for specifying output file path.
        *   ``.error`` -- flag for specifying error file path.
        *   ``.workdir`` -- flag for specifying path to working directory.
        *   ``.commands.submit`` -- submit command.
        *   ``.commands.check`` -- queue status check command.
        *   ``.commands.getid`` -- function extracting submitted job's ID from output of submit command.
        *   ``.commands.finished`` -- function checking if submitted job is finished. It should take a single string (job's ID) and return boolean.
        *   ``.commands.special.`` -- branch storing definitions of special |run| keyword arguments.

    See :meth:`call` for more technical details and examples.

    The *sleepstep* parameter defines how often the job is checked for being finished. It should be an integer value telling how many seconds should the interval between two checks last. If ``None``, the global default from ``config.sleepstep`` is copied.

    .. note::
        Usually job schedulers are configured in such a way that output of your job is captured somewhere else and copied to the location indicated by output flag when the job is finished. Because of that it is not possible to have a peek at your output while your job is running (for example to see if your calculation is going well). This limitation can be worked around with ``[Job].settings.runscript.stdout_redirect``. If set to ``True``, the output redirection will not be handled by a job scheduler, but built in the runscript using shell redirection ``>``. That forces the output file to be created directly in *workdir* and updated live as the job proceeds.
    """
    def __init__(self, grid='auto', sleepstep=None, **kwargs):
        JobRunner.__init__(self, **kwargs)
        self.sleepstep = sleepstep or config.sleepstep

        if grid == 'auto':
            self.settings = self._autodetect()
        elif grid in config.gridrunner:
            self.settings = config.gridrunner[grid]
            with open(os.devnull, 'wb') as null:
                try:
                    subprocess.call([self.settings.commands.submit, '--version'], stdout=null, stderr=null)
                except OSError:
                    raise PlamsError('GridRunner: %s command not found' %
                    self.settings.commands.submit)
        elif isinstance(grid, Settings):
            self.settings = grid
        else:
            raise PlamsError("GridRunner: invalid 'grid' argument. 'grid' should be either a Settings instance (see documentations for details) or a string occurring in config.gridrunner or 'auto' for autodetection")



    def call(self, runscript, workdir, out, err, runflags, **kwargs):
        """call(runscript, workdir, out, err, runflags, **kwargs)
        Submit *runscript* to the job scheduler with *workdir* as the working directory. Redirect output and error streams to *out* and *err*, respectively. *runflags* stores submit command options.

        The submit command has the following structure. Underscores denote spaces, parts in pointy brackets correspond to ``settings`` entries, parts in curly brackets to :meth:`call` arguments, square brackets contain optional parts::

            <.commands.submit>_<.workdir>_{workdir}_<.error>_{err}[_<.output>_{out}][FLAGS]_{runscript}

        Output part is added if *out* is not ``None``. This is handled automatically based on ``.runscript.stdout_redirect`` value in job's ``settings``.

        ``FLAGS`` part is built based on *runflags* argument, which is a dictionary storing |run| keyword arguments. For every *(key,value)* pair in *runflags* the string ``_-key_value`` is appended **unless** *key* is a special key occurring in ``.commands.special.``. In that case ``_<.commands.special.key>value`` is used (mind the lack of space in between!). For example, a |Settings| instance defining interaction with SLURM job scheduler stored in ``config.gridrunner.slurm`` has the following entries::

            .workdir = '-D'
            .output  = '-o'
            .error   = '-e'
            .special.nodes    = '-N '
            .special.walltime = '-t '
            .special.queue    = '-p '
            .commands.submit  = 'sbatch'
            .commands.check   = 'squeue -j '

        The submit command produced in such case::

            >>> gr = GridRunner(parallel=True, maxjobs=4, grid='slurm')
            >>> j.run(jobrunner=gr, queue='short', nodes=2, J='something', O='')

        will be::

            sbatch -D {workdir} -e {err} -o {out} -p short -N 2 -J something -O  {runscript}

        In some job schedulers some flags don't have a short form with semantics ``-key value``. For example, in SLURM the flag ``--nodefile=value`` have a short form ``-F value``, but the flag ``--export=value`` does not. One can still use such a flag using special keys mechanism::

            >>> gr = GridRunner(parallel=True, maxjobs=4, grid='slurm')
            >>> gr.settings.special.export = '--export='
            >>> j.run(jobrunner=gr, queue='short', export='value')
            sbatch -D {workdir} -e {err} -o {out} -p short --export=value {runscript}

        The submit command produced in the way explained above is then executed and returned output is used to determine submitted job's ID. Function stored in ``.commands.getid`` is used for that purpose, it should take one string (whole output) and return a string with job's ID.

        Now the method waits for the job to finish. Every ``sleepstep`` seconds it queries the job scheduler using following algorithm:
            *   if a key ``finished`` exists in ``.commands.`` it is used. It should be a function taking job's ID and returning ``True`` or ``False``.
            *   otherwise a string stored in ``.commands.check`` is concatenated with job's ID (with no space between) and such command is executed. Nonzero exit status indicates that job is no longer in job scheduler hence it is finished.

        Since it is difficult (on some systems even impossible) to automatically obtain job's exit code, the returned value is always 0. From |run| perspective it means that a job executed with |GridRunner| is never *crashed*.

        .. note::
            This method is used automatically during |run| and should never be explicitly called in your script.
        """

        s = self.settings
        cmd = ' '.join([s.commands.submit, s.workdir, workdir, s.error, err])
        if out is not None:
            cmd += ' '+s.output+' '+out
        for k,v in runflags.items():
            if k in s.special:
                cmd += ' '+s.special[k]+str(v)
            else:
                cmd += ' -'+k+' '+str(v)
        cmd += ' ' + opj(workdir,runscript)

        log('Submitting %s with command %s' % (runscript, cmd), 5)
        subout = subprocess.check_output(cmd.split(' '))
        subout = string(subout)
        log('Output of submit command: %s' % subout, 5)
        jobid = s.commands.getid(subout)
        log('%s submitted successfully as job %s' % (runscript, jobid), 3)

        if 'finished' in s.commands:
            while not s.commands.finished(jobid):
                time.sleep(self.sleepstep)
        else:
            checkcmd = (s.commands.check + jobid).split(' ')
            retcode = 0
            while retcode == 0:
                time.sleep(self.sleepstep)
                with open(os.devnull, 'wb') as null:
                    retcode = subprocess.call(checkcmd, stdout=null, stderr=null)

        log('Execution of %s finished' % runscript, 5)
        return 0



    def _autodetect(self):
        """Try to autodetect the type of job scheduler.

        The autodetection mechanism is very simple. For each entry in ``config.gridrunner`` the submit command followed by ``--version`` is executed (for example ``qsub --version``). If the execution was successful (which is indicated by exit code 0) the job scheduler of corresponding type is present on the system and it is chosen. So if there are multiple different job schedulers installed, only one is picked -- the one which "name" (indicated by a key in ``config.gridrunner``) is first in the lexicographical order.

        Returned value is one of ``config.gridrunner`` branches. If autodetection was not successful, an exception is raised.
        """
        for grid in config.gridrunner:
            with open(os.devnull, 'wb') as null:
                try:
                    retcode = subprocess.call([config.gridrunner[grid].commands.submit, '--version'], stdout=null, stderr=null)
                except OSError: continue
                if retcode == 0:
                    log('Grid type autodetected as ' + grid, 5)
                    return config.gridrunner[grid]
        raise PlamsError('GridRunner: Failed to autodetect grid type')
