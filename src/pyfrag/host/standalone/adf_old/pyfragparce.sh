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

Mark=(scf xc basis BECKEGRID EPRINT OCCUPATIONS ZLMFIT)
SMark=(numericalquality relativistic charge symmetry unrestricted PRINT STOFIT)


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
        echo "$item."${optionarray[0]}"="${optionarray[@]:1}
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
      echo $sitem"="${soptionarray[@]:1}
    else
      echo $sitem"=True"
    fi
  fi
done
}

input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

sed -n '/^JOBSUB$/,/^JOBSUB END$/{//!p;}' $input > jobsub.txt
sed -n '/^ADF$/,/^ADF END$/{//!p;}' $input > adf.txt
sed -n '/^PyFrag$/,/^PyFrag END$/{//!p;}' $input > pyfrag.txt
sed -n '/^fragment1 EXTRA$/,/^fragment1 EXTRA END$/{//!p;}' $input > fragment1_EXTRA.txt
sed -n '/^fragment2 EXTRA$/,/^fragment2 EXTRA END$/{//!p;}' $input > fragment2_EXTRA.txt
sed -n '/^complex EXTRA$/,/^complex EXTRA END$/{//!p;}' $input > complex_EXTRA.txt


submit="source \$HOSTPYFRAG/src/pyfrag/venv_runner.sh && run_pyfrag_cluster \$HOSTPYFRAG/src/pyfrag/cluster_runner.py --backend adf_old --job-dir \$PWD --input-file adfinputfile"


jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
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
