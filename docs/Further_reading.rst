Further Reading
===============

Whole Time Monitor
------------------

User can monitor the whole calculation process by using this command: ::

   pyfrag -x consist job.in

In this mode, every a period of time (default set is 20 seconds), a new data will be collected and updated in the form of webpage. In the meantime, if the original input is changed, a window will pop up to ask if you want to resummit job. If you agree to restart the job, a new input file will be submitted and started again.

In this case, you need specify the time interval to check the new result in /pyfrag/.pyfragrc  ::

   export JOBCHECK="20"

and in the /pyfrag/util/configure.py ::

   RESULTCHECK="20"
