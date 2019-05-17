Installation
============

Activation Strain Analysis (ASA) Module of PyFrag 2019
------------------------------------------------------
The user may choose to only install the part of the program needed to perform the Activation Strain Analysis (ASA) based on Activation Strain Model (ASM). Note that Python3 is needed to run this program. The ASA can be performed using a variety of quantum chemical software packages, including:  ADF_, Gaussian_, Orca_ and Turbomole_, given a series of coordinate from the potential energy surface is provided.

To install the ASA module of PyFrag 2019, the user must complete the following step. Go to your host machine (supercomputer or cluster), open a terminal and run the following command:

``curl -L -o install_alone.sh https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_alone.sh``

``bash install_alone.sh``

To run a simple test, open a terminal window on your host machine, make a directory, enter into that directory and run the following command to download the job input file (job.in) and coordinate file (molecule.xyz):

``curl -L -o job.in``
``https://raw.githubusercontent.com/sunxb05/PyFrag/master/host/standalone/adf_new/example/job.in``

``curl -L -o molecule.xyz``
``https://raw.githubusercontent.com/sunxb05/PyFrag/master/host/standalone/adf_new/example/molecule.xyz``

Change the ircpath and the submit information, such as the number of nodes and wall time, located in job.in using vim or any other text editor according to your situation, and run:

``pyfrag job.in``

The user can also download the module for either ADF, Gaussian, Orca, and Turbomole separately from PyFrag standalone_ and run it as a normal python code in your laptop or desktop. An input sample is provided in the example folder along with the source code file.


The Complete PyFrag 2019 Package
--------------------------------
The entire PyFrag 2019 package is only compatible with ADF at the moment. . For optimal use of PyFrag 2019, one part of the program is installed on the users’ local machine and the second part is installed on the users’ host machine (supercomputer or cluster) where the heavy computational jobs is running. The user must ensure to transport their public key to your host machine to allow the communication between your local and host machine. The following installation bash script (install_local.sh, install_host.sh) is was designed to make the installation process as simple as possible. However, for the advanced user, if a different configuration of the program is desired, please read the detailed comments in the installation bash script and set up the program accordingly.
To install and test PyFrag 2019, the user must perform the following three steps:


1) Go to your local machine (your laptop or desktop), open a terminal window and run the following command on your terminal:

``xcode-select --install``

``curl -L -o install_local.sh https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_local.sh``

``bash install_local.sh``


2) Go to your host machine (supercomputer or cluster), open a terminal window and run the following command:

``curl -L -o install_host.sh https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_host.sh``

``bash install_host.sh``

3)  Open a terminal window on your local machine, make a directory, enter into that directory and run the following command:

``curl -L -o job.in https://raw.githubusercontent.com/sunxb05/PyFrag/master/example/job.in``

Change the submit information, such as the number of nodes and wall time, located in job.in using vim or any other text editor, and run:

``pyfrag job.in``

To obtain the latest information about your job, the user can run:

``pyfrag -x summary job.in``


.. _PyFrag 2008: http://www.few.vu.nl/~xsn800/Home.html
.. _standalone: https://github.com/sunxb05/PyFrag/tree/master/host/standalone
.. _PyFrag 2019: https://sunxb05.github.io/pyfrag/
.. _Gaussian:   http://gaussian.com
.. _ADF:       https://www.scm.com
.. _Orca:      http://www.orcahome.de/orcanews.htm
.. _Turbomole: http://www.turbomole.com
