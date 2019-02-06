# title: PyFrag


## 1. Project summary

The PyFrag program is specially designed to facilitates the study of reaction mechanism in a more efficient and user-friendly way. It is an expansion of a popular program also named by PyFrag in our group. More information can be found in "https://sunxb05.github.io/pyfrag/". It automates the process of finding transition states, potential energy surface by using one simple input file. It follows by an activation strain analysis on the energy profile to characterize the feature of the reaction mechanism and gain insights into the overall reaction energies. Moreover, users can have an real-time monitoring of the running process via a webpage which vividly displays the updated data in the form of videos and figures and, if necessary, user can rerun the job immediately from where it stops. In this way, the three respects of computational chemistry–job management, data management and analysis management can all be contained in a framework and thus allow chemists to focus on the interpretation and creation work rather than waste time and energy on the finding and processing of massive data.

## 2. Installation

### 1 Prerequisite


1. ADF  ()     # computational engine, more information can be found in "https://www.scm.com". Noted the set-up should be handled properly so that adfinput can be called in terminal.
2. qmflow      # workflow manage engine, you can install qmflow following "https://github.com/SCM-NV/qmflows", however, it is needed to be replaced by the modified version which is located in ```/data/qmflows```
3. bokeh       # to show the website, more information can be found in "https://bokeh.pydata.org/en/latest/".
4. fabric      # to connect remote machine and transfer files, can be installed by pip ```pip install Fabric==1.12.2```


### 2 PyFrag

Download the resource of pyfrag which comprises two part, the local files and the remote files. The remote files ```/host``` include all files in the directory of host, which should be put in the remote machine where the heavy computation actually happens.

### 3 SETUP:

The path of bin of pyfrag should be put in the .bashrc or .profile in order to run pyfrag anywhere you want, something like:


export PYFRAGHOME="/Users/xiaobo/gitpyfrag"

export PATH=$PYFRAGHOME/bin:$PATH


Similarly, these related infomation should also be put in the .bashrc or .profile in your remote machine.

export HOSTPYFRAG='/home/x2sun/bin/host'

export QMWORKS='~/miniconda3/envs/qmworks'

export USERNAME='x2sun'



The basic setup is located in .pyfragrc, including the directory of videos of geometried generated in the optimization process, the local server, which means you need to set up your local websever service, and so on.

export PYFRAGVIDEO="/Users/xiaobo/Sites/video"

export PYFRAGHOST="http://localhost/~xiaobo/video"

export JOBCHECK="20"

export REMOTEBASE="/home/x2sun/pyfragtest_2"

export RESULTCHECK="20"

export HOSTPYFRAG='/home/x2sun/bin/host'


The configure for fabfile is located at utils directory, such as:

USERNAME = 'x2sun'

HOSTNAME = 'cartesius.surfsara.nl'

RESULTCHECK="20"



## 3 Usage

User can type pyfrag -h to see all the commands that can be used in this program, which will show something like:
Usage: /Users/xiaobo/gitpyfrag/bin/pyfrag [-h] [-s] [-x command]  [...]

       -h          : print this information
       -s          : run job quietly
       -x command  : start the executable named command
                   : command include restart, which restart job
                   : end, which terminate job
                   : check, which check the latest jobs information
                   : restart, which restart a job after it is stoped
                   : summary, which summarize all job result after jobs finished
                   : default command is pyfrag itself
The example command is like as follow, in which job.in is job input
/Users/xiaobo/gitpyfrag/bin/pyfrag job.in
or
/Users/xiaobo/gitpyfrag/bin/pyfrag -x restart job.in
or
/Users/xiaobo/gitpyfrag/bin/pyfrag -s -x summary job.in

## 4. Notice

### 1 Each job should be given a unique name and unique directory, because the result will be stored in the new directory named by the job name.

python3 -m pip install --user  --no-cache-dir  git+https://github.com/sunxb05/PyFrag@master#egg=qmworks-0.0.1

