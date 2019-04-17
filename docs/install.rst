Installation
============

Activation Strain Analysis
--------------------------
User could choose to install only part of the program to perform an activation strain analysis based on Activation Strain Model(ASM), which also known as the distortion/interaction model. The old PyFrag program (PyFrag2008_) can perform similar analysis using ADF as computational engine. This analysis now can be carried out on ADF_, Gaussian_, Orca_ and Turbomole_, given a series of coordinate from the potential energy surface is available.
To install the activation strain analysis module of PyFrag2019, user need do the following step:

Go to your host machine (supercomputer or cluster), open a terminal and run the following command on your terminal.

``curl -L -o install_alone.sh  https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_alone.sh``
``bash install_alone.sh``

User can also download the module for either ADF, Gaussian, Orca and Turbomole seperatedly from PyFrag standalone_ and run it as a normal python code in your laptop or desktop.
An input sample is in the example folder along with source code file.


Whole Package
-------------
The whole package only work for ADF right now. One part of the program is installed in your local machine, and the second part is installed in your host machine (suppercomputer or cluster) where the heavy computational jobs is running. So make sure you transport your public key to your host machine to allow the communication between your local and host machine. The following intallation bash script (install_local.sh, install_host.sh) is supposed to make the installation process as simple as possible. However, for the advanced user, if you want a different configuration of the program, please read the detailed comments in the intallation bash script and set up accordingly.

To install and test pyfrag, user need do three steps:

1) Go to your local machine (your laptop or desktop), open a terminal and run the following command on your terminal.

``curl -L -o install_local.sh  https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_local.sh``
``bash install_local.sh``


2) Go to your host machine (supercomputer or cluster), open a terminal and run the following command on your terminal.

``curl -L -o install_host.sh  https://raw.githubusercontent.com/sunxb05/PyFrag/master/install_host.sh``
``bash install_host.sh``

3) Open a terminal in your local machine, make a directory, enter into that directory and run the following command.

``curl -L -o job.in  https://raw.githubusercontent.com/sunxb05/PyFrag/master/example/job.in``

Change the submit information like the number of node and wall time in job.in using vim or any editor, and run:

``pyfrag job.in``

if your want to know the latest information about your job, just run:

``pyfrag -x summary job.in``


_PyFrag2008: http://www.few.vu.nl/~xsn800/Home.html
_standalone: https://github.com/sunxb05/PyFrag/tree/master/host/standalone
_PyFrag2019: https://sunxb05.github.io/pyfrag/
_Gaussian:   http://gaussian.com
_ADF:       https://www.scm.com
_Orca:      http://www.orcahome.de/orcanews.htm
_Turbomole: http://www.turbomole.com
