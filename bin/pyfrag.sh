#!/bin/sh

source $PYFRAGHOME/.pyfragbashrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"


#set up result dir
if  [ -e $JOBDIR/result ]; then
  rm $JOBDIR/result
fi
cp $PYFRAGHOME/data/result  $JOBDIR


#set for video directory for local server
if  [ -e $VIDEODIR ]; then
  rm $VIDEODIR
fi
mkdir $VIDEODIR
cp $PYFRAGHOME/data/video/index.html  $VIDEODIR


#save all information about job for future check


fab -f $PYFRAGHOME/utils/jobsub/fabfile.py deploy:$JOBDIR,$JOBNAME,$REMOTEDIR &
fab -f $PYFRAGHOME/utils/chresult/fabfile.py deploy:$REMOTEDIR/result,$JOBDIR/result,$REMOTEDIR,$JOBSTATE &


Mark="True"
while [ $Mark = "True" ]; do
  m1=$(md5sum "$JOBPATH")
  sleep $JOBCHECK
  m2=$(md5sum "$JOBPATH")
  if [ "$m1" != "$m2" ] ; then
    osascript $PYFRAGHOME/utils/buttons.scpt $JOBDIR $REMOTEDIR
  fi
  Mark=$(<$JOBSTATE)
done
