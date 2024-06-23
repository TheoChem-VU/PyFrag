Installation
============

Dependencies: in order to run pyfrag, the following programs should be installed:

   1. ADF_         (computational engine)
   2. qmworks_     (workflow management engine)
   3. bokeh_       (web server)
   4. fabric_      (data transfer program)


**Notes:**

- ADF_ is a popular DFT compuational program and its instalation is pretty staightford. However, after installation, the set-up of ADF_ should be handled properly so that adfinput can be called in terminal.
- qmworks_ is an extentison of workflow management python library qmflows_ developed by  F. Zapata, L. Ridder and B. F. van Beek which provide more functions than qmworks and more information can be find in the qmflows document_. The installation of qmworks_ is a bit tricky and is elaberated later.
- bokeh_ is an interactive visualization library that targets modern web browsers for presentation. It can be installed by ``pip install bokeh``
- fabric_ is a high level Python library designed to execute shell commands remotely over SSH, yielding useful Python objects in return. We use old version of fabirc which can be installed by ``pip install Fabric==1.12.2``.



.. _ADF: https://www.scm.com
.. _qmworks: https://github.com/sunxb05/PyFrag/tree/master/qmworks
.. _qmflows: https://github.com/SCM-NV/qmflows
.. _document: https://qmflows.readthedocs.io/en/latest/
.. _bokeh: https://bokeh.pydata.org/en/latest/
.. _fabric: http://www.fabfile.org




For a more detailed description of **instalation** read the documentation

.. toctree::
   local_machine
   remote_machine
