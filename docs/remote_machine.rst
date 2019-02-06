Remote machine
==============


Installation qmworks
--------------------

  - conda installation. Type in your console the following command:
    ```bash
     wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    ```

  - Then add miniconda to your path
    ```bash
     ./miniconda.sh -b -p $HOME/miniconda
    ```

  - Create new virtual environment
    ```bash
     conda create -q -n qmflows python=3.5
    ```

  - Install dependecies
    ```bash
     conda install --name qmflows -c anaconda hdf5
     conda install --name qmflows -c https://conda.anaconda.org/rdkit rdkit
    ```

  - Start environment
    ```bash
     source activate qmflows
    ```

  - install **qmflows** dependencies
    ```bash
     pip install https://github.com/SCM-NV/qmflows/tarball/master#egg=qmflows https://github.com/SCM-NV/plams/tarball/master#egg=plams --upgrade
    ```
### You are ready to start!



Similarly, these related infomation should also be put in the .bashrc or .profile in your remote machine.

    export HOSTPYFRAG='/home/x2sun/bin/host'
    export QMWORKS='~/miniconda3/envs/qmworks'
    export USERNAME='x2sun'
