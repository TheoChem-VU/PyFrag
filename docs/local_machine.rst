Local machine
=============

Download the resource of pyfrag_ and put it in a special directory like `bin` and then do the following set up.


The path of bin of pyfrag should be put in the .bashrc or .profile in your local machine in order to run pyfrag anywhere you want, something like: ::

  export PYFRAGHOME="/Users/yourname/PyFrag-master"
  (where you put the source code)
  export PATH=$PYFRAGHOME/bin:$PATH

The basic setup is located in .pyfragrc in the source code directory PyFrag-master, where you can set up the directory of videos of geometried generated in the optimization process, the local server, which means you need to set up your local websever service, and so on. For example ::

  export PYFRAGVIDEO="/Users/yourname/Sites/video"
  (where the video of geometrical optimization will be saved, be sure your local server is open)
  export PYFRAGHOST="http://localhost/video"
  (directory of video of your local website)
  export JOBCHECK="20"
  (time interval to check if your job input is changed)
  export REMOTEBASE="/home/pyfragtest_2"
  (the directory in your host machine where all jobs ins saved)
  export RESULTCHECK="20"
  (time interval to check if result is changed)
  export PYFRAGSOURCE='/home/bin/host'
  (where the code of host is stored in your host machine)


You need to modify the co configure for fabfile is located at utils directory, such as: ::

  USERNAME = 'x2sun'
  (the name of your account in your host machine)
  HOSTNAME = 'cartesius.surfsara.nl'
  (the address of your host machine)
  RESULTCHECK="20"
  (the same set as before)

.. _pyfrag: https://github.com/sunxb05/PyFrag
