#!/bin/sh

filename="/Users/xiaobo/Desktop/pyfrag/job.py"
jobstate="/Users/xiaobo/Desktop/pyfrag/result/jobstate.txt"

cp /Users/xiaobo/Dropbox/codes/bin/fab_subjob/jobstate.txt /Users/xiaobo/Desktop/pyfrag/result



fab -f $PYFRAGHOME/utils/jobsub/fabfile.py deploy:$JOBDIR,$JOBNAME,$REMOTEDIR &
fab -f $PYFRAGHOME/utils/chresult/fabfile.py deploy:$REMOTEDIR/result,$JOBDIR/result,$REMOTEDIR,$JOBSTATE &


Mark="True"
while [ $Mark = "True" ]; do
  m1=$(md5sum "$JOBPATH")
  # md5sum is computationally expensive, so check only once every 10 seconds
  sleep $JOBCHECK
  m2=$(md5sum "$JOBPATH")
  if [ "$m1" != "$m2" ] ; then
    osascript $PYFRAGHOME/utils/buttons.scpt $JOBDIR $REMOTEDIR
  fi
  Mark=$(<$JOBSTATE)
done
