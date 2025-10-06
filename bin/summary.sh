#!/bin/sh

kill $(ps aux | grep '[b]okeh' | awk '{print $2}') 2>/dev/null

source $PYFRAGHOME/.pyfragrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"



cp $JOBSTATE $JOBDIR/result


# Every a time interval, connect and check if new result is generated. If so, cp it to local machine.
/usr/local/bin/fab -f $PYFRAGHOME/utils/summary/fabfile.py deploy:$REMOTEDIR/result,$JOBDIR/result,$REMOTEDIR,$JOBSTATE,"${JOBNAME%.*}" &
