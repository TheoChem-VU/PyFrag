#!/bin/bash

#jobinfo:Submitted batch job 5354270
#runstate: JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) or
#          5354270     ...........................
#jobstate: True or False

function jobstate {
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"
process="$SCRIPTPATH/jobinfo.txt"
proID=`grep 'Submitted batch job' $process |awk '{print $4}'`
squeue -u $USERNAME > $SCRIPTPATH/result/runstate.txt
sleep 1
runstate="$SCRIPTPATH/result/runstate.txt"
jobstate=`grep  $proID  $runstate   |awk 'NF{ print $NF }'`
if   [ -z "$jobstate" ] ; then
printf 'False' > $SCRIPTPATH/result/jobstate.txt
fi
}


# function result_pro {
# cycle=$*
# eners=`grep -En 'GEOMETRY CONV|current energy' $cycle |awk '{print $3}'`
# enech=`grep 'energy change'            $cycle |awk '{print $3}'`
# gmaxs=`grep 'constrained gradient max' $cycle |awk '{print $4}'`
# grmss=`grep 'constrained gradient rms' $cycle |awk '{print $4}'`
# smaxs=`grep 'cart. step max'           $cycle |awk '{print $4}'`
# srmss=`grep 'cart. step rms'           $cycle |awk '{print $4}'`
# ener=`echo $eners |awk '{printf"%14.6f",627.509541*$1}'`
# printf "%15.4f %15.4f %15.6f %15.6f %15.6f %15.6f\n" $ener $enech $gmaxs $grmss $smaxs $srmss >  converge.txt
# grep -A 200 ' Coordinates' $cycle | grep -B 200 ' Number of elements' | grep -v ' ====' | grep -v '     Atom'| grep -v '                   X' | grep -v ' ------'| grep -v ' Number of elements' | grep -v ' Coordinates' > geometry.xyz
# }

function result_pro {
cycle=$*
eners=`grep -En 'GEOMETRY CONV|current energy' $cycle |awk '{print $3}'`
ener=`echo $eners |awk '{printf"%14.6f",-627.509541*$1}'`

enech_1=`grep 'energy change'            $cycle |awk '{print $3}'`
enech_2=`grep 'energy change'            $cycle |awk '{print $4}'`
enech_3=`grep 'energy change'            $cycle |awk '{print $5}'`

gmax_1=`grep 'constrained gradient max' $cycle |awk '{print $4}'`
gmax_2=`grep 'constrained gradient max' $cycle |awk '{print $5}'`
gmax_3=`grep 'constrained gradient max' $cycle |awk '{print $6}'`

grms_1=`grep 'constrained gradient rms' $cycle |awk '{print $4}'`
grms_2=`grep 'constrained gradient rms' $cycle |awk '{print $5}'`
grms_3=`grep 'constrained gradient rms' $cycle |awk '{print $6}'`

smax_1=`grep 'cart. step max'           $cycle |awk '{print $4}'`
smax_2=`grep 'cart. step max'           $cycle |awk '{print $5}'`
smax_3=`grep 'cart. step max'           $cycle |awk '{print $6}'`

srms_1=`grep 'cart. step rms'           $cycle |awk '{print $4}'`
srms_2=`grep 'cart. step rms'           $cycle |awk '{print $5}'`
srms_3=`grep 'cart. step rms'           $cycle |awk '{print $6}'`

# printf "%4.4f,%6.6f,%6.6f,%s,%6.6f,%6.6f,%s,%6.6f,%6.6f,%s,%6.6f,%6.6f,%s,%6.6f,%6.6f,%s\n" $ener $enech_1 $enech_2 $enech_3 $gmax_1 $gmax_2 $gmax_3 $grms_1 $grms_2 $grms_3 $smax_1 $smax_2 $smax_3 $srms_1 $srms_2 $srms_3>  converge.txt

printf "%4.4f %6.6f %6.6f %s %6.6f %6.6f %s %6.6f %6.6f %s %6.6f %6.6f %s %6.6f %6.6f %s\n" $ener $enech_1 $enech_2 $enech_3 $gmax_1 $gmax_2 $gmax_3 $grms_1 $grms_2 $grms_3 $smax_1 $smax_2 $smax_3 $srms_1 $srms_2 $srms_3>  converge.txt

grep -A 200 ' Coordinates' $cycle | grep -B 200 ' Number of elements' | grep -v ' ====' | grep -v '     Atom'| grep -v '                   X' | grep -v ' ------'| grep -v ' Number of elements' | grep -v ' Coordinates' > geometry.xyz
}




function result_final {
output=$*
step=`grep -En 'Geometry Convergence after Step  ' $output | cut -d":" -f 1`
steps=`echo $step |wc -w`
echo $steps
lstep='Geometry Convergence after Step   '$steps
grep -A 200 "$lstep" $output | grep -B 200 ' Number of elements of the density matrix on this node '    >   slot.txt
}



function compareconverge {
output=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

if [ ! -f $SCRIPTPATH/$2 ]; then
    cp  $1    $2
    cp  $1 ../../result/$3$1
elif cmp -s "$1" "$2"; then
    # printf " 0 "  >>   ../../result/result.txt
    # echo "0"
    true
else
    cp  $1    $2
    cp  $1 ../../result/$3$1
fi
}


function compare {
output=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

if [ ! -f $SCRIPTPATH/$2 ]; then
    cp  $1    $2
    cp  $1 ../../result/$3$1
    # printf ' 1 '  >>   ../../result/result.txt
elif cmp -s "$1" "$2"; then
    # echo "0"
    true
    # printf " 0 "  >>   ../../result/result.txt
else
    cp  $1    $2
    cp  $1 ../../result/$3$1
    # printf " 1 "   >>   ../../result/result.txt
fi
}


function irclog {
output=$*

tac $output | grep ' Coordinates in Geometry Cycle' -m 1 -B 9999 | tac | grep -B 100 'CORORT' | grep -v 'CORORT' | grep -v ' Coordinates'| grep -v 'Atom'    >   geometry.xyz
}




REMOTEDIR="$( cd "$(dirname "$1")" ; pwd -P )"
RESULTDIR=$REMOTEDIR/result

jobstate

# touch $RESULTDIR/result.txt
# echo 0 > $RESULTDIR/result.txt

unset -v latest
for file in "$REMOTEDIR"/plams*; do
  [[ $file -nt $latest ]] && latest=$file
done



if  [ -e $latest ]; then
  # if  [ -e $latest/*log ]; then
  #   cd $latest
  #   cp *log $RESULTDIR
  # fi
  for dname in r1 r2 rc ts p irc
  do
    if  [ -e $latest/$dname ]; then
        cd $latest/$dname
          if [ -e *out ]; then
            result_final *out
            result_pro slot.txt
            compareconverge  converge.txt _converge.txt $dname
            irclog logfile
            compare  geometry.xyz _geometry.xyz $dname
          fi
        cd $latest
    fi
  done
fi

if [ -e "$REMOTEDIR"/pyfrag* ]; then
  cp "$REMOTEDIR"/pyfrag* $RESULTDIR
  # printf " 1 "   >>   $RESULTDIR/result.txt
fi
