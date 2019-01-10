# -*- coding: utf-8 -*-

__author__ = 'xiaobo'


import os, re, sys
import time
from datetime import datetime
from fabric.api import *

# def arguments(path):
#   print (path)
#   return path

# _REMOTE_BASE_DIR = '/home/x2sun/pytest_1'
# _FILE            =  'job.py'
# current_path     =  '/Users/xiaobo/Desktop/pyfrag'
# current_path       =  sys.argv[1]
# exec("current_path" + "=" + '"' + sys.argv[1] + '"')

# def deploy():
#     with cd(_REMOTE_BASE_DIR):
#        run('bash /home/x2sun/bin/killjob.sh')
#        time.sleep(60)
#        run('rm -rf 1 2 3 4 5 *out adfinputfile jobinfo.txt job.py result sub')
#        with lcd(current_path):
#           local('bash /Users/xiaobo/Dropbox/codes/bin/fab_resub/restart.sh')

def deploy(path, argue):
      with lcd(path):
         print (argue)
         local('bash $PYFRAGHOME/utils/haha.sh')
# def deploy(path):
#       with lcd(path):
#          local('bash $PYFRAGHOME/utils/haha.sh')
