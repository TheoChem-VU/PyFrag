#!/bin/sh

source $PYFRAGHOME/.pyfragrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"



cp $JOBSTATE $JOBDIR/result

# connect to remote machine and send job input and start to parce and submit job.
fab -f $PYFRAGHOME/utils/jobsub/fabfile.py deploy:$JOBDIR,$JOBNAME,$REMOTEDIR &

# Every a time interval, connect and check if new result is generated. If so, cp it to local machine.
fab -f $PYFRAGHOME/utils/chresult/fabfile.py deploy:$REMOTEDIR/result,$JOBDIR/result,$REMOTEDIR,$JOBSTATE,"${JOBNAME%.*}" &


#Check if input file is changed, if changed, there are three options: resubmit, end computation or do nothing.
Mark="True"
while [ $Mark = "True" ]; do
  m1=$(md5sum "$JOBPATH")
  sleep $JOBCHECK
  m2=$(md5sum "$JOBPATH")
  if [ "$m1" != "$m2" ] ; then
    osascript $PYFRAGHOME/utils/buttons.scpt $JOBDIR $REMOTEDIR $JOBNAME
  fi
  Mark=$(<$JOBSTATE)
done
