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

Mark=(scf xc basis BeckeGrid EPRINT OCCUPATIONS ZLMFIT Properties REMOVEFRAGORBITALS)
SMark=(numericalquality relativistic symmetry unrestricted STOFIT Print REMOVEALLFRAGVIRTUALS)


for item in ${Mark[*]}
do
grep  -iA 20 "$item" $adf | grep -m 1 -iB 20 'end' | grep -iv 'end' | grep -iv "$item" > "$item".txt
  if [ -s "$item".txt ]; then
    argue="$item".txt
    while read -r line
    do
      if [ "$line" != "" ]; then
        option="$line"
        optionarray=( $option )
        echo "adf.""$item."${optionarray[0]}"="${optionarray[@]:1}
      fi
    done < "$argue"
  fi
rm "$item".txt
done

for sitem in ${SMark[*]}
do
  soption=`grep -iw "$sitem" $adf`
  if [ ! -z "$soption" ]; then
    soptionarray=( $soption )
    if [[ ! -z ${soptionarray[@]:1} ]]; then
      echo "adf."$sitem"="${soptionarray[@]:1}
    else
      echo "adf."$sitem"=True"
    fi
  fi
done



soption=`grep -iw "charge" $adf`
if [ ! -z "$soption" ]; then
  soptionarray=( $soption )
  if [[ ! -z ${soptionarray[@]:1} ]]; then
    echo "ams.system.charge="${soptionarray[@]:1}
  fi
fi

}



input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

grep  -A 200 'JOBSUB' $input | grep -B 200 'JOBSUB END' | grep -v 'JOBSUB' | grep -v 'JOBSUB END' > jobsub.txt
grep  -A 200 'ADF' $input | grep -B 200 'ADF END' | grep -v 'ADF' | grep -v 'ADF END' > adf.txt
grep  -A 200 'PyFrag' $input | grep -B 200 'PyFrag END' | grep -v 'PyFrag' | grep -v 'PyFrag END' > pyfrag.txt
grep  -A 200 'fragment1 EXTRA' $input | grep -B 200 'fragment1 EXTRA END' | grep -v 'fragment1 EXTRA' | grep -v 'fragment1 EXTRA END' > fragment1_EXTRA.txt
grep  -A 200 'fragment2 EXTRA' $input | grep -B 200 'fragment2 EXTRA END' | grep -v 'fragment2 EXTRA' | grep -v 'fragment2 EXTRA END' > fragment2_EXTRA.txt
grep  -A 200 'complex EXTRA' $input | grep -B 200 'complex EXTRA END' | grep -v 'complex EXTRA' | grep -v 'complex EXTRA END' > complex_EXTRA.txt


submit="amspython \$HOSTPYFRAG/standalone/adf_new/PyFrag.py \\"
subadfinputfile="--adfinputfile "$SCRIPTPATH/"adfinputfile \\"

echo -n ""                                                   > ./sub
jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub
echo $subadfinputfile                                       >> ./sub
echo -n ""                                                   > ./adfinputfile
adfargue  adf.txt                                           >> ./adfinputfile

extraOption=(fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt)

for extraItem in ${extraOption[*]}
do
  if [ -s $extraItem ]; then
    echo -n ""                                               > ./${extraItem%.txt}
    adfargue $extraItem                                     >> ./${extraItem%.txt}
    subItem="--"${extraItem%.txt}" $SCRIPTPATH/${extraItem%.txt} \\"
    echo $subItem                                           >> ./sub
  fi
done

rm jobsub.txt adf.txt pyfrag.txt fragment1_EXTRA.txt fragment2_EXTRA.txt complex_EXTRA.txt
