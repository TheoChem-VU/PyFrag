Installation
============

Download the resource of pyfrag which comprises two part, the local files and the remote files. The remote files ```/host``` include all files in the directory of host, which should be put in the remote machine where the heavy computation actually happens.


Dependencies: in order to run pyfrag, the following programs should be installed: ::
   1. ADF_         (computational engine)
   2. qmflows_     (workflow manage engine, however, it is needed to be replaced by the modified version)
   3. bokeh_       (to show the website.)
   4. fabric_      (to connect remote machine and transfer files.)


**Notes:**

- Old version of Fabirc can be installed by ``pip install Fabric==1.12.2``.
- Noted the set-up of ADF_ should be handled properly so that adfinput can be called in terminal.



.. _ADF: https://www.scm.com
.. _qmflows: https://github.com/SCM-NV/qmflows
.. _bokeh: https://bokeh.pydata.org/en/latest/
.. _fabric: http://www.fabfile.org




For a more detailed description of **instalation** read the documentation

.. toctree::
   local_machine
   remote_machine
