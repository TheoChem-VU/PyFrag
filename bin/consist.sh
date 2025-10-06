#!/bin/sh


source $PYFRAGHOME/.pyfragrc

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOBPATH=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
REMOTEDIR=$REMOTEBASE/"${JOBNAME%.*}"


#set up result dir
#Each job should be given a unique name, because the result will be stored in the new directory named by the job name.
random=$(( $RANDOM % 20000 + 5000 ))
if  [ -e $JOBDIR/result ]; then
  rm -rf $JOBDIR/result
fi
# Including the file of jobstate.txt, the content of which is "True", means jobs is running
cp -r $PYFRAGHOME/data/result  $JOBDIR
echo $random > $JOBDIR/result/stocks/port.txt


#set for video directory for local server, named by the job name(without appenddix)
if  [ -e $VIDEODIR ]; then
  rm -r $VIDEODIR
fi
mkdir $VIDEODIR
cp  -r $PYFRAGHOME/data/video/csv2html $VIDEODIR



#save all information about job for future check

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

kill -9 $(ps -p $PPID -o ppid=)
