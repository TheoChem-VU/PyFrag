#!/bin/sh

source $PYFRAGHOME/.pyfragrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"




/usr/local/bin/fab -f $PYFRAGHOME/utils/termjob/fabfile.py deploy:$JOBDIR,$REMOTEDIR
/usr/local/bin/fab -f $PYFRAGHOME/utils/resub_simple/fabfile.py deploy:$JOBDIR,$REMOTEDIR,$JOBNAME &
