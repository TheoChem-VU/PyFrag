# -*- coding: utf-8 -*-

__author__ = 'xiaobo'

'''
Deployment toolkit.
'''
import time
import os
from os import sys, path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configure import *
from fabric.api import *

env.user = USERNAME
env.hosts = [HOSTNAME]



def deploy(REMOTERESULT,CURRENTRESULT,REMOTEDIR,JOBSTATE,JOBNAME):
    Mark="True"
    while Mark=="True":
      time.sleep(RESULTCHECK)
      with cd(REMOTEDIR):
        run('bash $HOSTPYFRAG/result/result.sh')
        get('%s/%s' % (REMOTERESULT, "jobstate.txt"), CURRENTRESULT)
        r=int(run('bash $HOSTPYFRAG/result/change.sh ./result/result.txt'))
        if r >=1:
#          print ('retrun argue', r)
          get('%s/%s' % (REMOTERESULT, "*"), CURRENTRESULT)
          with lcd(CURRENTRESULT):
            local('bash $PYFRAGHOME/process/figure.sh')
            local('bash $PYFRAGHOME/server/web.sh %s' % JOBNAME)
            local('bash $PYFRAGHOME/server/bokeh.sh')
      with open(str(JOBSTATE), "r") as fd:
        all_words = fd.read().splitlines()
      fd.close()
      Mark=all_words[0]
