Further Reading
===============

Whole Time Monitor
------------------

User can have full time monitor the whole calculation process: ::

``pyfrag -x consist job.in``

In this mode, every a period of time, a new data will be collected and update in the form of webpage.
In the meantime, if the original input is changed, a window will pop up to ask if you want to resummit job.
If you agree to restart job, a new input file will be submitted and start again.

In this case, you need to do more configuration in the .pyfragrc in the source code directory of PyFrag and the util/Pyrag, where you need to specify the time interval to check the new result. The default set is 20 second. For example ::

  export JOBCHECK="20"
  (time interval to check if your job input is changed)


  RESULTCHECK="20"
  (the same set as before)
