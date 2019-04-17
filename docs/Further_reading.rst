Further Reading
===============

Whole Time Monitor
------------------

User can monitor the whole calculation process by using this command: ::

   pyfrag -x consist job.in

In this mode, every a period of time, a new data will be collected and updated in the form of webpage. In the meantime, if the original input is changed, a window will pop up to ask if you want to resummit job. If you agree to restart the job, a new input file will be submitted and started again.

In this case, you need do more configuration source code directory /pyfrag/.pyfragrc and in the /pyfrag/util/configure.py, where you need to specify the time interval to check the new result. The default set is 20 seconds. For example ::
  #in /pyfrag/.pyfragrc
  export JOBCHECK="20"

  #in /pyfrag/util/configure.py
  RESULTCHECK="20"

