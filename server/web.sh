#!/bin/sh


function webedit {

htmlfile=$RESULTDIR/stocks/templates/index.html

for dir in r1geometry r2geometry rcgeometry tsgeometry pgeometry ircgeometry
do
    if [ -d $RESULTDIR/$dir ]; then
      # dirName="$(basename "$dir")"
      cp -r    $RESULTDIR/$dir       $PYFRAGVIDEO/$JOBNAME
      cp $PYFRAGHOME/data/video/video_index.html  $PYFRAGVIDEO/$JOBNAME/$dir/index.html
      echo '    <div class="videoWrapper">'                                                                                                    >> $htmlfile
      echo '    <iframe width="820" height="615" src="'$PYFRAGHOST'/'$JOBNAME'/'$dir'/" frameborder="0" allowfullscreen></iframe>'             >> $htmlfile
      echo '    </div>'                                                                                                                        >> $htmlfile
      echo '    <p><strong>'$dir':</strong> Video of updated geometrical structure in the calculation.</p>'                                    >> $htmlfile
    fi
done

echo '    <iframe width="820" height="615" src="'$PYFRAGHOST'/'$JOBNAME'/csv2html/" frameborder="0" allowfullscreen></iframe>'>> $htmlfile
echo '</div>'                                                                                                                                  >> $htmlfile
echo '{% endblock %}'                                                                                                                          >> $htmlfile
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

function pyfragfig {
if  [ -e $RESULTDIR/pyfragdefault.txt ]; then
  mv $RESULTDIR/pyfragdefault.txt   $RESULTDIR/stocks/daily/pyfrag.txt
  cd $RESULTDIR/stocks
  python3 change.py $RESULTDIR/stocks/daily/pyfrag.txt
fi
}

function geotable {
if  [ -e $RESULTDIR ]; then
  for dname in r1geometry.xyz r2geometry.xyz rcgeometry.xyz tsgeometry.xyz pgeometry.xyz ircgeometry.xyz
  do
    fileName="${dname%geometry.xyz}"
    if  [ -f $RESULTDIR/$dname ]; then
          if [ -s $RESULTDIR/$dname ]; then
            cd $RESULTDIR
            python3 $PYFRAGHOME/server/geotable.py $RESULTDIR/$dname $RESULTDIR/"$fileName"converge.txt $PYFRAGVIDEO/$JOBNAME/csv2html/data/coor.csv
          fi
    fi
  done
fi
}

function pyfragtable {

if  [ -e $RESULTDIR/pyfrag.txt ]; then
  python3 $PYFRAGHOME/server/pyfragtable.py $RESULTDIR/pyfrag.txt $PYFRAGVIDEO/$JOBNAME/csv2html/data/coor.csv
fi
}


RESULTDIR="$( pwd -P )"
JOBNAME=$*
# JOBDIR="$( cd "$(dirname "$RESULTDIR")" ; pwd -P )"
# JOBNAME="$(basename "$JOBDIR")"
# initiate index.html for results
cp $PYFRAGHOME/data/video/pre_index.html   $RESULTDIR/stocks/templates/index.html
# edit index.html to update result, including title and figures
webedit
#figure for process monitor
convergefig

cp $PYFRAGHOME/data/video/csv2html/data/coor.csv $PYFRAGVIDEO/$JOBNAME/csv2html/data/coor.csv
#table for latest coordinate
geotable

#figure for pyfrag
pyfragfig

#table for pyfrag
pyfragtable

# if  [ -e $RESULTDIR/pyfrag* ]; then
#   cp $RESULTDIR/pyfrag*   $RESULTDIR/stocks/daily/pyfrag.txt
#   cd $RESULTDIR/stocks
#   python3 change.py $RESULTDIR/stocks/daily/pyfrag.txt
# fi
