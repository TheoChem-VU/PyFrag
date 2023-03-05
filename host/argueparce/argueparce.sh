function jobsubargue {
jobsub="$1"
while read -r line
do
    if [ "$line" != "" ]; then
      name="$line"
      echo "$name"
    fi
done < "$jobsub"
}

function pyfragargue {
pyfrag="$1"
while read -r line
do
    if [ "$line" != "" ]; then
      name="$line"
      echo "--$name \\"
    fi
done < "$pyfrag"
}

function adfargue {
adf=$*

Mark=(scf xc basis BECKEGRID)
SMark=(numericalquality relativistic charge symmetry Restart)


for item in ${Mark[*]}
do
  option=`grep  -iA 20 "$item" $adf | grep -m 1 -iB 20 'end' | grep -iv 'end' | grep -iv "$item"`
  optionarray=( $option )
  i=0
  j=1
  while test "$j" -le "${#optionarray[@]}"
  do
    echo $item"."${optionarray[$i]}"="${optionarray[$j]}
    i=`expr $i + 2`
    j=`expr $j + 2`
  done
done

for sitem in ${SMark[*]}
do
  soption=`grep -i "$sitem" $adf`
  if [ ! -z "$soption" ]; then
    soptionarray=( $soption )
    echo $sitem"="${soptionarray[@]:1}
  fi
done

grep  -iA 20 "EPRINT" $adf | grep -m 1 -iB 20 'END' | grep -iv 'END' | grep -iv "EPRINT"  > eprint.txt


if [ -s eprint.txt ]; then
  eprint=eprint.txt
  while read -r line
  do
      if [ "$line" != "" ]; then
        option="$line"
        optionarray=( $option )
        echo "eprint."${optionarray[0]}"="${optionarray[@]:1}
      fi
  done < "$eprint"
fi

rm eprint.txt


grep  -iA 20 "tsrc" $adf | grep -m 1 -iB 20 'end' | grep -iv 'end' | grep -iv "tsrc"  > tsrc.txt


if [ -s tsrc.txt ]; then
  tsrc=tsrc.txt
  while read -r line
  do
      if [ "$line" != "" ]; then
        option="$line"
        optionarray=( $option )
        echo "tsrc."${optionarray[0]}"="${optionarray[@]:1}
      fi
  done < "$tsrc"
fi

rm tsrc.txt


grep  -iA 20 "geometry" $adf | grep -m 1 -iB 20 'end' | grep -iv 'end' | grep -iv "geometry"  > geometry.txt


if [ -s geometry.txt ]; then
  geometry=geometry.txt
  while read -r line
  do
      if [ "$line" != "" ]; then
        option="$line"
        optionarray=( $option )
        length=${#optionarray[@]}
        if [ ${length} == 3 ]; then
          echo "geometry."${optionarray[0]}"."${optionarray[1]}"="${optionarray[@]:2}
        else
          echo "geometry."${optionarray[0]}"="${optionarray[@]:1}
        fi
      fi
  done < "$geometry"
fi

rm geometry.txt
}


input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

sed -n '/^JOBSUB$/,/^JOBSUB END$/{//!p;}' $input > jobsub.txt
sed -n '/^ADF$/,/^ADF END$/{//!p;}' $input > adf.txt
sed -n '/^PyFrag$/,/^PyFrag END$/{//!p;}' $input > pyfrag.txt
sed -n '/^Geometrycoor$/,/^Geometrycoor END$/{//!p;}' $input > coor.xyz
sed -n '/^R1 EXTRA$/,/^R1 EXTRA END$/{//!p;}' $input > R1_EXTRA.txt
sed -n '/^R2 EXTRA$/,/^R2 EXTRA END$/{//!p;}' $input > R2_EXTRA.txt
sed -n '/^RC EXTRA$/,/^RC EXTRA END$/{//!p;}' $input > RC_EXTRA.txt
sed -n '/^TS EXTRA$/,/^TS EXTRA END$/{//!p;}' $input > TS_EXTRA.txt
sed -n '/^P EXTRA$/,/^P EXTRA END$/{//!p;}' $input > P_EXTRA.txt
sed -n '/^IR EXTRA$/,/^IR EXTRA END$/{//!p;}' $input > IRC_EXTRA.txt
sed -n '/^IR_1 EXTRA$/,/^IR_1 EXTRA END$/{//!p;}' $input > IRC_1_EXTRA.txt
sed -n '/^fragment1 EXTRA$/,/^fragment1 EXTRA END$/{//!p;}' $input > fragment1_EXTRA.txt
sed -n '/^fragment2 EXTRA$/,/^fragment2 EXTRA END$/{//!p;}' $input > fragment2_EXTRA.txt
sed -n '/^complex EXTRA$/,/^complex EXTRA END$/{//!p;}' $input > complex_EXTRA.txt


submit="$QMWORKS/bin/python3 $HOSTPYFRAG/job.py \\"
subadfinputfile="--adfinputfile "$SCRIPTPATH/"adfinputfile \\"


jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
python3 $HOSTPYFRAG/argueparce/coorargue.py coor.xyz        >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub
echo $subadfinputfile                                       >> ./sub
adfargue  adf.txt                                           >> ./adfinputfile

extraOption=(R1_EXTRA.txt R2_EXTRA.txt RC_EXTRA.txt TS_EXTRA.txt P_EXTRA.txt IRC_EXTRA.txt IRC_1_EXTRA.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt)

for extraItem in ${extraOption[*]}
do
  if [ -s $extraItem ]; then
    adfargue $extraItem                                     >> ./${extraItem%.txt}
    subItem="--"${extraItem%.txt}" $SCRIPTPATH/${extraItem%.txt} \\"
    echo $subItem                                           >> ./sub
  fi
done

rm jobsub.txt adf.txt pyfrag.txt coor.xyz R1_EXTRA.txt R2_EXTRA.txt RC_EXTRA.txt TS_EXTRA.txt P_EXTRA.txt IRC_EXTRA.txt IRC_1_EXTRA.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt
mkdir result
touch ./result/rcgeometry.xyz
cp $HOSTPYFRAG/argueparce/jobstate.txt  ./result
