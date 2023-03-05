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

input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

sed -n '/^PyFrag$/,/^PyFrag END$/{//!p;}' $input > pyfrag.txt

submit="python3 \$HOSTPYFRAG/standalone/adf_openorb/PyFrag.py \\"

echo $submit                                                >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub

rm  pyfrag.txt
