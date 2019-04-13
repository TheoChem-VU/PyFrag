#!/usr/local/bin/python3
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

env.user  = USERNAME
env.hosts = [HOSTNAME]


def deploy(REMOTERESULT,CURRENTRESULT,REMOTEDIR,JOBSTATE,JOBNAME):
  with cd(REMOTEDIR):
    run('bash $HOSTPYFRAG/result/result.sh')
    get('%s/%s' % (REMOTERESULT, "*"), CURRENTRESULT)
    with lcd(CURRENTRESULT):
      local('bash $PYFRAGHOME/process/figure.sh')
      local('bash $PYFRAGHOME/server/web.sh %s' % JOBNAME)
      local('bash $PYFRAGHOME/utils/summary/bokeh.sh')
