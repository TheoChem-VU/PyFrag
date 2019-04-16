Installation
============

To install and test pyfrag, user need do three steps:

1) Go to your local machine (your laptop or desktop), open a terminal and run the following command on your terminal.

``curl -L -o install_local.sh  https://github.com/sunxb05/PyFrag/install_local.sh``
``bash install_local.sh``


2) Go to your host machine (supercomputer or cluster), open a terminal and run the following command on your terminal.

``curl -L -o install_host.sh  https://github.com/sunxb05/PyFrag/install_host.sh``
``bash install_host.sh``

3) Open a terminal in your local machine, make a directory, enter into that directory and run the following command.

``curl -L -o job.in  https://github.com/sunxb05/PyFrag/blob/master/example/job.in``

change the submit information in job.in using vim or any editor, and run:

``pyfrag job.in``

if your want to know the latest information about your job, just run:

``pyfrag -x summary job.in``
