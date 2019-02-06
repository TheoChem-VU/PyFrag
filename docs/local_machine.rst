Local machine
=============


Download the resource of pyfrag which comprises two part, the local files and the remote files. The remote files ```/host``` include all files in the directory of host, which should be put in the remote machine where the heavy computation actually happens.


SETUP:
------

The path of bin of pyfrag should be put in the .bashrc or .profile in order to run pyfrag anywhere you want, something like: ::


  export PYFRAGHOME="/Users/xiaobo/gitpyfrag"
  export PATH=$PYFRAGHOME/bin:$PATH


The basic setup is located in .pyfragrc, including the directory of videos of geometried generated in the optimization process, the local server, which means you need to set up your local websever service, and so on. For example ::

  export PYFRAGVIDEO="/Users/xiaobo/Sites/video"
  export PYFRAGHOST="http://localhost/~xiaobo/video"
  export JOBCHECK="20"
  export REMOTEBASE="/home/x2sun/pyfragtest_2"
  export RESULTCHECK="20"
  export HOSTPYFRAG='/home/x2sun/bin/host'


The configure for fabfile is located at utils directory, such as: ::

  USERNAME = 'x2sun'
  HOSTNAME = 'cartesius.surfsara.nl'
  RESULTCHECK="20"
