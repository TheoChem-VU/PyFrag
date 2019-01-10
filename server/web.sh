#!/bin/sh

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"
JOBNAME=`basename "$1"`
JOB=$JOBDIR/$JOBNAME
VIDEODIR=$PYFRAGVIDEO/"${JOBNAME%.*}"
JOBSTATE=$JOBDIR/result/jobstate.txt


function webedit {
loginput=~/Desktop/pyfrag/result/*log
htmlfile=/Users/xiaobo/Sites/node/bokeh/examples/app/stocks/templates/index.html
logfile=~/Desktop/pyfrag/result/log.txt

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
for dir in ~/Desktop/pyfrag/result/*
do
    if [ -d $dir ]; then
      dirName="$(basename "$dir")"
      cp -r $dir /Users/xiaobo/Sites/video
      cp /Users/xiaobo/Dropbox/codes/bin/fab_chresult/index.html /Users/xiaobo/Sites/video/$dirName
      echo '    <div class="videoWrapper">'    >> $htmlfile
      echo '    <iframe width="820" height="615" src="http://localhost/~xiaobo/video/'$dirName'/" frameborder="0" allowfullscreen></iframe>'             >> $htmlfile
      echo '    </div>'                        >> $htmlfile
    fi
done

echo '</div>'                                  >> $htmlfile
echo '{% endblock %}'                          >> $htmlfile

}



SCRIPTPATH="/Users/xiaobo/Desktop/pyfrag"


function convergefig {

if  [ -e $SCRIPTPATH/result ]; then
  for dname in r1converge.txt r2converge.txt rcconverge.txt tsconverge.txt pconverge.txt
  do
    fileName="${dname%converge.txt}"
    if  [ -f $SCRIPTPATH/result/$dname ]; then
          if [ -s $SCRIPTPATH/result/$dname ]; then
            cd $SCRIPTPATH/result
            conFile=$(<$dname)
            echo $conFile >> /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/daily/$fileName
            python3 /Users/xiaobo/Dropbox/codes/bin/fab_chresult/converge.py /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/daily/$fileName
          fi
    fi
  done
fi

}


SCRIPTPATH="/Users/xiaobo/Desktop/pyfrag"

# initiate index.html for results
cp /Users/xiaobo/Dropbox/codes/bin/fab_chresult/bokeh/index.html   /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/templates


# edit index.html to update result, including title and figures
webedit

#figure for process monitor
convergefig


#figure for pyfrag

if  [ -e $SCRIPTPATH/result/pyfrag* ]; then
  cp ~/Desktop/pyfrag/result/pyfrag*   /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/daily/pyfrag.txt
  python3 /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/change.py /Users/xiaobo/Sites/node/bokeh/examples/app/stocks/daily/pyfrag.txt
fi
