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

Mark=(scf xc basis)
SMark=(numericalquality relativistic charge symmetry)


for item in ${Mark[*]}
do
  option=`grep  -A 20 "$item" $adf | grep -m 1 -B 20 'end' | grep -v 'end' | grep -v "$item"`
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
  soption=`grep "$sitem" $adf`
  soptionarray=( $soption )
  echo $sitem"="${soptionarray[@]:1}
done
}


input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

grep  -A 200 'JOBSUB' $input | grep -B 200 'JOBSUB END' | grep -v 'JOBSUB' | grep -v 'JOBSUB END' > jobsub.txt

grep  -A 200 'ADF' $input | grep -B 200 'ADF END' | grep -v 'ADF' | grep -v 'ADF END' > adf.txt

grep  -A 200 'PyFrag' $input | grep -B 200 'PyFrag END' | grep -v 'PyFrag' | grep -v 'PyFrag END' > pyfrag.txt

grep  -A 200 'Geometrycoor' $input | grep -B 200 'Geometrycoor END' | grep -v 'Geometrycoor' | grep -v 'Geometrycoor END' > coor.xyz

#jobname=`grep '#SBATCH -J' $input |awk '{print $1}'`
submit="~/miniconda3/envs/qmworks/bin/python3 /home/x2sun/bin/job.py \\"
subadfinputfile="--adfinputfile "$SCRIPTPATH/"adfinputfile"


jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
python3 /home/x2sun/bin/argueparce/coorargue.py coor.xyz    >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub
echo $subadfinputfile                                       >> ./sub

adfargue  adf.txt                                           >> ./adfinputfile

rm jobsub.txt adf.txt pyfrag.txt coor.xyz
rm -r result
mkdir result
touch ./result/rcgeometry.xyz
cp /home/x2sun/bin/argueparce/jobstate.txt  ./result
