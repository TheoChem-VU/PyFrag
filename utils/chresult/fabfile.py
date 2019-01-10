# -*- coding: utf-8 -*-

__author__ = 'xiaobo'

'''
Deployment toolkit.
'''

import os, re
import time
from datetime import datetime
from fabric.api import *
from ..configure import *




def deploy(REMOTERESULT,CURRENTRESULT,REMOTEDIR,JOBSTATE):
    Mark="True"
    while Mark=="True":
      time.sleep(RESULTCHECK)
      with cd(REMOTEDIR):
        run('bash $HOSTPYFRAG/result/result.sh')
        get('%s/%s' % (REMOTERESULT, "jobstate.txt"), CURRENTRESULT)
        r=int(run('bash $HOSTPYFRAG/result/change.sh result.txt'))
        if r >=1:
#          print ('retrun argue', r)
          get('%s/%s' % (REMOTERESULT, "*"), CURRENTRESULT)
          with lcd(CURRENTRESULT):
            local('bash $PYFRAGHOME/process/figure.sh')
            local('bash $PYFRAGHOME/server/web.sh')
            local('bash $PYFRAGHOME/server/bokeh.sh')
      with open(str(JOBSTATE), "r") as fd:
        all_words = fd.read().splitlines()
      fd.close()
      Mark=all_words[0]
