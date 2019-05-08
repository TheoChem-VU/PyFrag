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
SMark=(numericalquality relativistic charge symmetry)


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
}


input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

grep  -A 200 'JOBSUB' $input | grep -B 200 'JOBSUB END' | grep -v 'JOBSUB' | grep -v 'JOBSUB END' > jobsub.txt
grep  -A 200 'ADF' $input | grep -B 200 'ADF END' | grep -v 'ADF' | grep -v 'ADF END' > adf.txt
grep  -A 200 'PyFrag' $input | grep -B 200 'PyFrag END' | grep -v 'PyFrag' | grep -v 'PyFrag END' > pyfrag.txt
grep  -A 200 'fragment1 EXTRA' $input | grep -B 200 'fragment1 EXTRA END' | grep -v 'fragment1 EXTRA' | grep -v 'fragment1 EXTRA END' > fragment1_EXTRA.txt
grep  -A 200 'fragment2 EXTRA' $input | grep -B 200 'fragment2 EXTRA END' | grep -v 'fragment2 EXTRA' | grep -v 'fragment2 EXTRA END' > fragment2_EXTRA.txt
grep  -A 200 'complex EXTRA' $input | grep -B 200 'complex EXTRA END' | grep -v 'complex EXTRA' | grep -v 'complex EXTRA END' > complex_EXTRA.txt


submit="python3 \$HOSTPYFRAG/standalone/adf_new/PyFrag.py \\"
subadfinputfile="--adfinputfile "$SCRIPTPATH/"adfinputfile \\"


jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub
echo $subadfinputfile                                       >> ./sub
adfargue  adf.txt                                           >> ./adfinputfile

extraOption=(fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt)

for extraItem in ${extraOption[*]}
do
  if [ -s $extraItem ]; then
    adfargue $extraItem                                     >> ./${extraItem%.txt}
    subItem="--"${extraItem%.txt}" $SCRIPTPATH/${extraItem%.txt} \\"
    echo $subItem                                           >> ./sub
  fi
done

rm jobsub.txt adf.txt pyfrag.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt
