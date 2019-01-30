# -*- coding: utf-8 -*-

__author__ = 'xiaobo'

import os
from os import sys, path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configure import *
from fabric.api import *

env.user = USERNAME
env.hosts = [HOSTNAME]

def deploy(CURRENT_DIR, FILE, REMOTE_DIR):

    # run('mkdir %s' % REMOTE_DIR)
    run('if [ ! -d %s ]; then mkdir %s; fi' % (REMOTE_DIR, REMOTE_DIR))
    put('%s/%s' % (CURRENT_DIR, FILE), REMOTE_DIR)

    with cd(REMOTE_DIR):
       run('bash $HOSTPYFRAG/argueparce/argueparce.sh %s' % FILE)
       run('source activate qmworks')
       run('sbatch sub > jobinfo.txt')
