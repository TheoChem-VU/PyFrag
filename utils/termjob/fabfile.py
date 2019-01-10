# -*- coding: utf-8 -*-

__author__ = 'xiaobo'


import os, re
from fabric.api import *
from ..configure import *



def deploy(REMOTEDIR):
    with cd(REMOTEDIR):
        run('bash $HOSTPYFRAG/result/killjob.sh')


