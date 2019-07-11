from __future__ import unicode_literals

import glob
import os
import shutil
import threading
import time
import types

from os.path import join as opj
from six import PY3
if PY3:
    import builtins
else:
    import __builtin__ as builtins


from .errors import PlamsError
from .settings import Settings

__all__ = ['init', 'finish', 'log', 'load', 'load_all', 'add_to_class', 'add_to_instance']


#===================================================================================================

def init(path=None, folder=None):
    """Initialize PLAMS environment. Create global ``config`` and default |JobManager|.

    An empty |Settings| instance is created and added to :mod:`public<__builtin__>` namespace as ``config``. Then it is populated with default settings by executing ``plams_defaults.py``. The following locations are used to search for the defaults file, in order of precedence:
        *   If ``$PLAMSDEFAULTS`` variable is in your environment and it points to a file, this file is used (executed as Python script).
        *   If ``$PLAMSHOME`` variable is in your environment and ``$PLAMSHOME/utils/plams_defaults.py`` exists, it is used.
        *   If ``$ADFHOME`` variable is in your environment and ``$ADFHOME/scripting/plams/utils/plams_defaults.py`` exists, it is used.
        *   Otherwise, the path ``../../../utils/plams_defaults.py`` relative to the current file (``common.py``) is checked. If defaults file is not found there, an exception is raised.

    Next, a |JobManager| instance is created as ``config.jm`` using *path* and *folder* to determine the main working directory. Settings used by this instance are directly linked from ``config.jobmanager``. If *path* is not supplied, the current directory is used. If *folder* is not supplied, the string ``plams.`` followed by PID of the current process is used.

    .. warning::
      This function **must** be called before any other PLAMS command can be executed. Trying to do anything without it results in a crash. See also :ref:`master-script`.
    """

    builtins.config = Settings()

    from os.path import isfile, expandvars, dirname
    if 'PLAMSDEFAULTS' in os.environ and isfile(expandvars('$PLAMSDEFAULTS')):
        defaults = expandvars('$PLAMSDEFAULTS')
    elif 'PLAMSHOME' in os.environ and isfile(opj(expandvars('$PLAMSHOME'), 'src', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$PLAMSHOME'), 'src', 'scm', 'plams', 'plams_defaults')
    elif 'ADFHOME' in os.environ and isfile(opj(expandvars('$ADFHOME'), 'scripting', 'plams', 'src', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$ADFHOME'), 'scripting', 'plams', 'src', 'scm', 'plams', 'plams_defaults')
    else:
        defaults = opj(dirname(dirname(__file__)), 'plams_defaults')
        if not isfile(defaults):
            raise PlamsError('plams_defaults not found, please set PLAMSDEFAULTS or PLAMSHOME in your environment')
    exec(compile(open(defaults).read(), defaults, 'exec'))


    from .jobmanager import JobManager
    config.jm = JobManager(config.jobmanager, path, folder)

    log('PLAMS running with Python %i' % (3 if PY3 else 2), 5)
    log('PLAMS environment initialized', 5)
    log('PLAMS working folder: %s' % config.jm.workdir, 1)


#===================================================================================================

def finish(otherJM=None):
    """Wait for all threads to finish and clean the environment.

    This function must be called at the end of your script for :ref:`cleaning` to take place. See :ref:`master-script` for details.

    If for some reason you use other job managers than the default one, they need to passed as *otherJM* list.
    """
    for thread in threading.enumerate():
        if thread.name == 'plamsthread':
            thread.join()

    config.jm._clean()
    if otherJM:
        for jm in otherJM:
            jm._clean()
    log('PLAMS environment cleaned up', 5)

    if config.erase_workdir is True:
        shutil.rmtree(config.jm.workdir)


#===================================================================================================

def load(filename):
    """Load previously saved job from ``.dill`` file. This is just a shortcut for |load_job| method of the default |JobManager| ``config.jm``."""
    return config.jm.load_job(filename)


#===================================================================================================

def load_all(path, jobmanager=None):
    """Load all jobs from *path*.

    This function works as a multiple execution of |load_job|. It searches for ``.dill`` files inside the directory given by *path*, yet not directly in it, but one level deeper. In other words, all files matching ``path/*/*.dill`` are used. That way a path to the main working folder of previously run script can be used to import all jobs run by that script.

    The purpose of this function is to provide quick and easy way of restarting a script that previously failed. Loading all successful jobs from the previous run prevents double work and allows the script to proceed directly to the place where it failed.

    Jobs are loaded using default job manager stored in ``config.jm``. If you wish to use a different one you can pass it as *jobmanager* argument of this function.

    Returned value is a dictionary containing all loaded jobs as values and absolute paths to ``.dill`` files as keys.
    """
    jm = jobmanager or config.jm
    loaded_jobs = {}
    for f in glob.glob(opj(path, '*', '*.dill')):
        loaded_jobs[f] = jm.load_job(f)
    return loaded_jobs


#===================================================================================================

_stdlock = threading.Lock()
_filelock = threading.Lock()

def log(message, level=0):
    """Log *message* with verbosity *level*.

    Logs are printed independently to both text file and standard output. If *level* is equal or lower than verbosity (defined by ``config.log.file`` or ``config.log.stdout``) the message is printed. Date and/or time can be added based on ``config.log.date`` and ``config.log.time``. All logging activity is thread safe.
    """
    if 'config' in vars(builtins):
        if level <= config.log.file or level <= config.log.stdout:
            message = str(message)
            prefix = ''
            if config.log.date:
                prefix += '%d.%m|'
            if config.log.time:
                prefix += '%H:%M:%S'
            if prefix:
                prefix = '[' + prefix.rstrip('|') + '] '
                message = time.strftime(prefix) + message
            if level <= config.log.stdout:
                with _stdlock:
                    print(message)
            if level <= config.log.file and 'jm' in config:
                with _filelock, open(config.jm.logfile, 'a') as f:
                    f.write(message + '\n')


#===================================================================================================

def add_to_class(classname):
    """Add decorated function as a method to the whole class *classname*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        @add_to_class(ADFResults)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

    After executing the above code all instances of ``ADFResults`` (even the ones created earlier) are enriched with ``get_energy`` method that can be invoked by::

        someadfresults.get_energy()

    The added method is visible from subclasses of *classname* so ``@add_to_class(Results)`` will also work in the above example.

    If *classname* is |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see :ref:`parallel`).
    """
    from .results import _restrict, _MetaResults
    def decorator(func):
        if isinstance(classname, _MetaResults):
            func = _restrict(func)
        setattr(classname, func.__name__, func)
    return decorator

#===================================================================================================

def add_to_instance(instance):
    """Add decorated function as a method to one particular *instance*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        results = myjob.run()

        @add_to_instance(results)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

        results.get_energy()

    The added method is visible only for one particular instance and it overrides any methods defined on class level or added with :func:`add_to_class` decorator.

    If *instance* is an instance of |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see :ref:`parallel`).
    """
    from .results import _restrict, Results
    def decorator(func):
        if isinstance(instance, Results):
            func = _restrict(func)
        func = types.MethodType(func, instance)
        setattr(instance, func.__func__.__name__, func)
    return decorator


#===================================================================================================

#remove me and all my calls after moving to Python3!!!
def string(s):
    if PY3 and isinstance(s, bytes):
        return s.decode()
    return s

