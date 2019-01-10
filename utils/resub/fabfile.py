# -*- coding: utf-8 -*-

__author__ = 'xiaobo'


import os, re, sys
import time
from fabric.api import *
from ..configure import *




def deploy(JOBDIR, REMOTEDIR):
    with cd(REMOTEDIR):
       run('bash $HOSTPYFRAG/result/killjob.sh')
       time.sleep(ESULTCHECK)
       run('rm -rf 1 2 3 4 5 *out adfinputfile jobinfo.txt job.py result sub')
       with lcd(JOBDIR):
          local('bash $PYFRAGHOME/bin/resubpy')
