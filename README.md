---
title: Pyfrag
date: 2019-01-07 17:47:22
---


## 1. PyFrag in combination with other programs

The PyFrag program is specially designed to facilitates the study of reaction mechanism in a more efficient and user-friendly way. It automates the process of finding transition states, potential energy surface by using one simple input file. It follows by an activation strain analysis on the energy profile to characterize the feature of the reaction mechanism and gain insights into the overall reaction energies. Moreover, users can have an real-time monitoring of the running process via a webpage which vividly displays the updated data in the form of videos and figures and, if necessary, user can rerun the job immediately from where it stops. In this way, the three respects of computational chemistry–job management, data management and analysis management can all be contained in a framework and thus allow chemists to focus on the interpretation and creation work rather than waste time and energy on the finding and processing of massive data.

## 2. PyFrag installation

# 1
Download the resource of pyfrag which comprise two part, the local files and the remote files. The remote files include all files in the directory of host, which should be put in the remote machine where the heavy computation actually happens.

# 2

SETUP:
The path of bin of pyfrag should be put in the .bashrc or .profile in order to run pyfrag anywhere you want, something like:
# Setting PATH for PyFrag
export PYFRAGHOME="/Users/xiaobo/gitpyfrag"
export PATH=$PYFRAGHOME/bin:$PATH
# Setting PATH for PyFrag

Similarly, these related infomation should also be put in the .bashrc or .profile in your remote machine.

export HOSTPYFRAG='/home/x2sun/bin/host'
export QMWORKS='~/miniconda3/envs/qmworks'



The basic setup is located in .pyfragrc, including the directory of videos of geometried generated in the optimization process, the local server, which means you need to set up your local websever service, and so on.

export PYFRAGVIDEO="/Users/xiaobo/Sites/video"
export PYFRAGHOST="http://localhost/~xiaobo/video"
export JOBCHECK="20"       #time interval set to check if input file is changed
export REMOTEBASE="/home/x2sun/pyfragtest"
export RESULTCHECK="20"       #time interval set to check if result is changed
export HOSTPYFRAG='/home/x2sun/bin/host'
export QMWORKS='~/miniconda3/envs/qmworks/bin/python3'
export USERNAME='x2sun'

The configure for fabfile is located at utils directory, such as:
USERNAME = 'x2sun'
HOSTNAME = 'cartesius.surfsara.nl'
HOSTPYFRAG='/home/x2sun/bin/host'
RESULTCHECK="20"    #time interval set to check if result is changed



As well, the same should be applied in your host machine, that is, the supercompuer or cluster where you run real compuatations. It should be something like:


Pyfrag can be splited into two part. One part should be installed in local machine, and the other should be put into the remote host machine (for example suppercommuter or cluster that do the heavy compuational jobs). This program relies on some of the following programs.

1. ADF  ()     # computational engine
2. qmflow      # workflow manage engine
3. bokeh       # to show the website
4. fabric      # to connect remote machine

## 5. Notice

# 1 Each job should be given a unique name, because the result will be stored in the new directory named by the job name.



