#Preview mode: no calculations are actually run, only inputs and runscripts are prepared.
config.preview = False

#Defines a 'quantum' of time used whenever some action needs to be repeated until a certain condition is met (i.e. querying queue manager about the job submitted to a grid system)
config.sleepstep = 5

#Allows to try to obtain results of failed/crashed jobs. Could be disasterous, so beware :)
config.ignore_failure = True

#If set to True, all threads started by JobRunner are daemon threads. Daemon threads are terminated when the main thread finishes and hence allow immediate end of the whole script when Ctrl-C is pressed.
config.daemon_threads = True

#if set to True, the entire main working folder is deleted at the end of script
config.erase_workdir = False

#======================= JobManager defaults ==========================================

#When two or more jobs have the same name they are renamed to [jobname].002. Defines number of digits of appended number
config.jobmanager.counter_len = 3

#Defines hashing method for testing if job was previously run. Currently supported values are 'input', 'runscript', 'input+runscript' and False/None
config.jobmanager.hashing = 'input'

#Removes all empty subdirectories in working directory at the end of the script
config.jobmanager.remove_empty_directories = True

#Defines what action to take when a particular jobfolder already exists in the filesystem. Possible values are: None (throw exception), 'remove', 'rename' (to *.old). Relevant only when running using an existing main working folder.
config.jobmanager.jobfolder_exists = None

#======================= Job defaults =================================================

#After job is finished, pickle job object to [jobname].dill
config.job.pickle = True

#Define which files produced by executed job should be kept. See manual for details and possible values
config.job.keep = 'all'
config.job.save = 'all'

#First line of produced runscripts
config.job.runscript.shebang = '#!/bin/sh'

#If set to True, '>[jobname].out' is used in runscript. Otherwise python takes care of output redirection (set it to True if you want to peek output of job submitted to grid)
config.job.runscript.stdout_redirect = False

#When files are imported into job's directory (by rerun prevention), they can be either copied or hardlinked, the following setting adjusts this behavior. On Windows it has no effect, files are always copied.
config.job.link_files = True

#======================= Log defaults =================================================

#Verbosity of log messages: 0:none  1:minimal  3:normal  5:verbose  7:extremely talkative
#Verbosity of log printed to .log file in working directory
config.log.file = 5
#Verbosity of log printed to standard output
config.log.stdout = 3
#Print time for each log event
config.log.time = True
#Print date for each log event
config.log.date = False

#======================= Default JobRunner ============================================

from .jobrunner import JobRunner
config.default_jobrunner = JobRunner(parallel=False)




#======================================================================================
#===================== GridRunner definitions==========================================
#======================================================================================

# GridRunner mechanism for testing if job is finished:
#if [...].commands.finished exists it is used to check if the job is finished. It should be a function that takes a single string (job_id) as an argument and returns True or False
#otherwise [...].commands.check is combined with job_id, executed as a subprocess and returned exit code is tested (nonzero return code indicates that job has finished)

def __getid(s):
    return ''.join([c for c in s if c.isdigit()])

config.gridrunner.pbs.workdir = '-d'
config.gridrunner.pbs.output  = '-o'
config.gridrunner.pbs.error   = '-e'
config.gridrunner.pbs.special.nodes    = '-l nodes='
config.gridrunner.pbs.special.walltime = '-l walltime='
config.gridrunner.pbs.special.queue    = '-q '
config.gridrunner.pbs.commands.submit  = 'qsub'
config.gridrunner.pbs.commands.check   = 'qstat '
config.gridrunner.pbs.commands.getid   = __getid

config.gridrunner.slurm.workdir = '-D'
config.gridrunner.slurm.output  = '-o'
config.gridrunner.slurm.error   = '-e'
config.gridrunner.slurm.special.nodes    = '-N '
config.gridrunner.slurm.special.walltime = '-t '
config.gridrunner.slurm.special.queue    = '-p '
config.gridrunner.slurm.commands.submit  = 'sbatch'
config.gridrunner.slurm.commands.check   = 'squeue -j '
config.gridrunner.slurm.commands.getid   = __getid

def __slurm_finished(jobid):
    try:
        import subprocess32 as subprocess
    except ImportError:
        import subprocess
    from .common import string
    import os
    cmd = ('squeue -j ' + jobid).split(' ')
    try:
        with open(os.devnull, 'wb') as null:
            out = subprocess.check_output(cmd, stderr=null)
    except subprocess.CalledProcessError:
        return True
    lines = string(out).split('\n')
    return not (len(lines) > 1 and lines[1].find(jobid) != -1)

config.gridrunner.slurm.commands.finished   = __slurm_finished

