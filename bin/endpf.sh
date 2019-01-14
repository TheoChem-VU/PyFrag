#!/bin/sh

filename="/Users/xiaobo/Desktop/pyfrag/job.py"
jobstate="/Users/xiaobo/Desktop/pyfrag/result/jobstate.txt"

cp /Users/xiaobo/Dropbox/codes/bin/fab_subjob/jobstate.txt /Users/xiaobo/Desktop/pyfrag/result

fab -f $PYFRAGHOME/utils/termjob/fabfile.py deploy:$REMOTEDIR &

