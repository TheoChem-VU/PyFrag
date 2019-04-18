Further Reading
===============

Whole Time Monitor
------------------

The user can monitor the entire calculation process by using this the following command: ::

   pyfrag -x consist job.in

In this mode, periodically (default set is 20 seconds), new data will be collected and updated in the form of webpage. In the meantime, if the original input is changed, a window will pop up to ask the user if they want to resubmit the job. If the user agrees to restart the job, a new input file will be submitted and started again and all other previous data will remain unchanged.

In this case, you need specify the time interval to check the new result in /pyfrag/.pyfragrc  ::

   export JOBCHECK="20"

and in the /pyfrag/util/configure.py ::

   RESULTCHECK="20"
