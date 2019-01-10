#!/bin/sh
JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOB=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt
RESULTDIR=$JOBDIR/result

function webedit {
loginput=$RESULTDIR/*log
htmlfile=$RESULTDIR/stocks/templates/index.html
logfile=$RESULTDIR/log.txt

touch $logfile
grep  'finished with status' $loginput  > $logfile

if [ -s $logfile ]; then
  loginfo=$logfile
else
  loginfo=$loginput
fi


while read line
do
echo -e  '<br>' $line                          >> $htmlfile
done < $loginfo


echo '    </p>'                                >> $htmlfile
echo '    {{ super() }}'                       >> $htmlfile
for dir in $RESULTDIR/*
do
    if [ -d $dir ]; then
      dirName="$(basename "$dir")"
      cp -r $dir $PYFRAGVIDEO
      cp $PYFRAGHOME/data/video/video_index.html  $PYFRAGVIDEO/$dirName/index.html
      echo '    <div class="videoWrapper">'    >> $htmlfile
      echo '    <iframe width="820" height="615" src="$PYFRAGHOST/'$dirName'/" frameborder="0" allowfullscreen></iframe>'             >> $htmlfile
      echo '    </div>'                        >> $htmlfile
    fi
done

echo '</div>'                                  >> $htmlfile
echo '{% endblock %}'                          >> $htmlfile
}


function convergefig {
if  [ -e $RESULTDIR ]; then
  for dname in r1converge.txt r2converge.txt rcconverge.txt tsconverge.txt pconverge.txt
  do
    fileName="${dname%converge.txt}"
    if  [ -f $RESULTDIR/$dname ]; then
          if [ -s $RESULTDIR/$dname ]; then
            cd $RESULTDIR
            conFile=$(<$dname)
            echo $conFile >> $RESULTDIR/stocks/daily/$fileName
            python3 $PYFRAGHOME/server/converge.py $RESULTDIR/stocks/daily/$fileName
          fi
    fi
  done
fi
}


# initiate index.html for results
cp $PYFRAGHOME/data/video/pre_index.html   $RESULTDIR/stocks/templates
# edit index.html to update result, including title and figures
webedit
#figure for process monitor
convergefig
#figure for pyfrag

if  [ -e $RESULTDIR/pyfrag* ]; then
  cp $RESULTDIR/pyfrag*   $RESULTDIR/stocks/daily/pyfrag.txt
  cd $RESULTDIR/stocks
  python3 change.py $RESULTDIR/stocks/daily/pyfrag.txt
fi
