process="jobinfo.txt"
proID=`grep 'Submitted batch job' $process |awk '{print $4}'`
scancel $proID

