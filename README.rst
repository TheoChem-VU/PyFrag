PyFrag
#######
See documentation_ for tutorials and documentation.

Motivation
==========
The PyFrag program is specially designed to facilitates the study of reaction mechanism in a more efficient and user-friendly way. It is an expansion of a popular program also named by PyFrag in our group. More information can be found in "https://sunxb05.github.io/pyfrag/". It automates the process of finding transition states, potential energy surface by using one simple input file. It follows by an activation strain analysis on the energy profile to characterize the feature of the reaction mechanism and gain insights into the overall reaction energies. Moreover, users can have an real-time monitoring of the running process via a webpage which vividly displays the updated data in the form of videos and figures and, if necessary, user can rerun the job immediately from where it stops. In this way, the three respects of computational chemistry–job management, data management and analysis management can all be contained in a framework and thus allow chemists to focus on the interpretation and creation work rather than waste time and energy on the finding and processing of massive data.

Description
===========
This library consists of a set of modules written in Python3 to
automate the following tasks:

 1. Input generation.
 2. Handle tasks dependencies (Noodles_).
 3. Advanced molecular manipulation capabilities with (rdkit_).
 4. Numerical data storage and manipulation (HDF5_).
 5. Jobs failure detection and recovery.
 6. Distribution in heterogeneous hardware platforms.

Tutorial and Examples
---------------------
A tutorial written as a jupyter-notebook is available from notebook_. You can
also access direclty more examples.


Installation
============


For installation, please read installation_.


.. _documentation: https://pyfragdocument.readthedocs.io/en/latest/
.. _installation: https://pyfragdocument.readthedocs.io/en/latest/install.html
.. _notebook: https://github.com/sunxb05/PyFrag/tree/master/jupyterNotebooks/
.. _installConda: http://conda.pydata.org/docs/install/quick.html
.. _Noodles: http://nlesc.github.io/noodles/
.. _HDF5: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _ADF: https://www.scm.com
.. _examples: https://github.com/SCM-NV/qmflows/tree/master/src/qmflows/examples
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _qd-example: https://github.com/SCM-NV/qmflows/blob/master/test/QD_input_examples

