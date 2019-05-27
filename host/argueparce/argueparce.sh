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
        echo "geometry."${optionarray[0]}"="${optionarray[@]:1}
      fi
  done < "$geometry"
fi

rm geometry.txt
}


input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

grep  -A 200 'JOBSUB' $input | grep -B 200 'JOBSUB END' | grep -v 'JOBSUB' | grep -v 'JOBSUB END' > jobsub.txt
grep  -A 200 'ADF' $input | grep -B 200 'ADF END' | grep -v 'ADF' | grep -v 'ADF END' > adf.txt
grep  -A 200 'PyFrag' $input | grep -B 200 'PyFrag END' | grep -v 'PyFrag' | grep -v 'PyFrag END' > pyfrag.txt
grep  -A 200 'Geometrycoor' $input | grep -B 200 'Geometrycoor END' | grep -v 'Geometrycoor' | grep -v 'Geometrycoor END' > coor.xyz
grep  -A 200 'R1 EXTRA' $input | grep -B 200 'R1 EXTRA END' | grep -v 'R1 EXTRA' | grep -v 'R1 EXTRA END' > R1_EXTRA.txt
grep  -A 200 'R2 EXTRA' $input | grep -B 200 'R2 EXTRA END' | grep -v 'R2 EXTRA' | grep -v 'R2 EXTRA END' > R2_EXTRA.txt
grep  -A 200 'RC EXTRA' $input | grep -B 200 'RC EXTRA END' | grep -v 'RC EXTRA' | grep -v 'RC EXTRA END' > RC_EXTRA.txt
grep  -A 200 'TS EXTRA' $input | grep -B 200 'TS EXTRA END' | grep -v 'TS EXTRA' | grep -v 'TS EXTRA END' > TS_EXTRA.txt
grep  -A 200 'P EXTRA' $input | grep -B 200 'P EXTRA END' | grep -v 'P EXTRA' | grep -v 'P EXTRA END'     > P_EXTRA.txt
grep  -A 200 'IR EXTRA' $input | grep -B 200 'IR EXTRA END' | grep -v 'IR EXTRA' | grep -v 'IR EXTRA END' > IRC_EXTRA.txt
grep  -A 200 'fragment1 EXTRA' $input | grep -B 200 'fragment1 EXTRA END' | grep -v 'fragment1 EXTRA' | grep -v 'fragment1 EXTRA END' > fragment1_EXTRA.txt
grep  -A 200 'fragment2 EXTRA' $input | grep -B 200 'fragment2 EXTRA END' | grep -v 'fragment2 EXTRA' | grep -v 'fragment2 EXTRA END' > fragment2_EXTRA.txt
grep  -A 200 'complex EXTRA' $input | grep -B 200 'complex EXTRA END' | grep -v 'complex EXTRA' | grep -v 'complex EXTRA END' > complex_EXTRA.txt


submit="$QMWORKS/bin/python3 $HOSTPYFRAG/job.py \\"
subadfinputfile="--adfinputfile "$SCRIPTPATH/"adfinputfile \\"


jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
python3 $HOSTPYFRAG/argueparce/coorargue.py coor.xyz        >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub
echo $subadfinputfile                                       >> ./sub
adfargue  adf.txt                                           >> ./adfinputfile

extraOption=(R1_EXTRA.txt R2_EXTRA.txt RC_EXTRA.txt TS_EXTRA.txt P_EXTRA.txt IRC_EXTRA.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt)

for extraItem in ${extraOption[*]}
do
  if [ -s $extraItem ]; then
    adfargue $extraItem                                     >> ./${extraItem%.txt}
    subItem="--"${extraItem%.txt}" $SCRIPTPATH/${extraItem%.txt} \\"
    echo $subItem                                           >> ./sub
  fi
done

rm jobsub.txt adf.txt pyfrag.txt coor.xyz R1_EXTRA.txt R2_EXTRA.txt RC_EXTRA.txt TS_EXTRA.txt P_EXTRA.txt IRC_EXTRA.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt
mkdir result
touch ./result/rcgeometry.xyz
cp $HOSTPYFRAG/argueparce/jobstate.txt  ./result
