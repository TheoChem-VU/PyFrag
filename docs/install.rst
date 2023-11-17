Installation
============
The user may choose to only install the part of the program needed to perform the Activation Strain Analysis (ASA) based on Activation Strain Model (ASM). Note that Python3 is needed to run this program. The ASA can be performed using a variety of quantum chemical software packages, including:  ADF_, Gaussian_, Orca_ and Turbomole_, given a series of coordinate from the potential energy surface is provided.

Activation Strain Analysis (ASA) Module of PyFrag 2019.0.1
----------------------------------------------------------
To install the latest version (development) of the ASA module of PyFrag 2019, go to your host machine (supercomputer or cluster), open a terminal and run the following command:

``curl -L -o install_alone.sh https://raw.githubusercontent.com/TheoChem-VU/PyFrag/development/install_alone.sh``

``bash install_alone.sh``

To run a simple test, open a terminal window on your host machine, make a directory, enter into that directory and run the following command to download the job input file (job.in) and coordinate file (molecule.xyz):

``curl -L -o new_ams_job.in``
``https://raw.githubusercontent.com/TheoChem-VU/PyFrag/development/host/standalone/adf_new/example/new_ams_job.in``

``curl -L -o molecule.xyz``
``https://raw.githubusercontent.com/TheoChem-VU/PyFrag/development/host/standalone/adf_new/example/molecule.xyz``

Change the ircpath and the submit information, such as the number of nodes and wall time, located in job.in using vim or any other text editor according to your situation, and run:

``pyfrag new_ams_job.in``

We recommend installing this version since it has the most recent bugfixes.

Activation Strain Analysis (ASA) Module of PyFrag 2019
------------------------------------------------------
To install the ASA module of PyFrag 2023, go to your host machine (supercomputer or cluster), open a terminal and run the following command:

``curl -L -o install_alone.sh https://raw.githubusercontent.com/TheoChem-VU/PyFrag/master/install_alone.sh``

``bash install_alone.sh``

To run a simple test, open a terminal window on your host machine, make a directory, enter into that directory and run the following command to download the job input file (job.in) and coordinate file (molecule.xyz):

``curl -L -o job.in``
``https://raw.githubusercontent.com/TheoChem-VU/PyFrag/master/host/standalone/adf_new/example/job.in``

``curl -L -o molecule.xyz``
``https://raw.githubusercontent.com/TheoChem-VU/PyFrag/master/host/standalone/adf_new/example/molecule.xyz``

Change the ircpath and the submit information, such as the number of nodes and wall time, located in job.in using vim or any other text editor according to your situation, and run:

``pyfrag job.in``


------------------------------------------------------

The user can also download the module for either ADF, Gaussian, Orca, and Turbomole separately from PyFrag standalone_ and run it as a normal python code in your laptop or desktop. An input sample is provided in the example folder along with the source code file.

.. _PyFrag 2008: http://www.few.vu.nl/~xsn800/Home.html
.. _standalone: https://github.com/TheoChem-VU/PyFrag/tree/master/host/standalone
.. _PyFrag 2019: https://sunxb05.github.io/pyfrag/
.. _Gaussian:   http://gaussian.com
.. _ADF:       https://www.scm.com
.. _Orca:      http://www.orcahome.de/orcanews.htm
.. _Turbomole: http://www.turbomole.com
