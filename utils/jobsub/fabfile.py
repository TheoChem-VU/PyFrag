# -*- coding: utf-8 -*-

__author__ = 'xiaobo'


import os, re
import time
from datetime import datetime
from fabric.api import *
from ..configure import *


def deploy(CURRENT_DIR, FILE, REMOTE_DIR):
    put('%s/%s' % (CURRENT_DIR, FILE), REMOTE_DIR)

    with cd(REMOTE_DIR):
       run('bash $HOSTPYFRAG/argueparce/argueparce.sh $FILE')
       run('source activate qmworks')
       run('sbatch sub > jobinfo.txt')
