process="jobinfo.txt"
RESULTDIR="$( pwd -P )"
if  [ -f $RESULTDIR/$process ]; then
  proID=`grep 'Submitted batch job' $process |awk '{print $4}'`
  if [ -z $proID ]; then
    echo "Job already finished"
  else
    scancel $proID
  fi
fi


