#!/bin/sh

source $PYFRAGHOME/.pyfragrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"



cp $JOBSTATE $JOBDIR/result

# connect to remote machine and send job input and start to parse and submit job.
/usr/local/bin/fab -f $PYFRAGHOME/utils/jobsub/fabfile.py deploy:$JOBDIR,$JOBNAME,$REMOTEDIR &
