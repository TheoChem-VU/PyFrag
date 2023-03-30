### General information ###
# This script is used create an parser input file for the PyFrag.py script
# The pyfrag input file (*not* the parser input file) is split into three main sections 
# and possibly extra sections depending of additional ADF/AMS settings for the two fragments and complex
#
# The three main sections are:
# JOBSUB contains sbatch information (including loading modules) and ends up in the "sub" file
# AMS    contains the ADF/AMS settings that is used for the two fragments and complex (common settings)
# PyFrag contains parser commands specifically for the PyFrag.py script and PyFrag related options
# 
# Additional sections are:
# fragment1 EXTRA contains ADF/AMS settings that is used for fragment 1 only 
# fragment2 EXTRA contains ADF/AMS settings that is used for fragment 2 only
# complex EXTRA   contains ADF/AMS settings that is used for the complex only
#
# This is only compatible with the AMS2020 and later versions of AMS! 

function jobsubargue {
# This function translates the jobsub.txt file into a string of arguments that ends up in the "sub" file
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
# This function translates the pyfrag.txt file into a string of arguments (format: "--parser_key") that ends up in the "sub" file
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

sed -n '/^JOBSUB$/,/^JOBSUB END$/{//!p;}' $input > jobsub.txt
sed -n '/^AMS$/,/^AMS END$/{//!p;}' $input > adfinputfile
sed -n '/^ADF$/,/^ADF END$/{//!p;}' $input > old_adfinputfile
sed -n '/^PyFrag$/,/^PyFrag END$/{//!p;}' $input > pyfrag.txt
sed -n '/^fragment1 EXTRA$/,/^fragment1 EXTRA END$/{//!p;}' $input > fragment1_EXTRA
sed -n '/^fragment2 EXTRA$/,/^fragment2 EXTRA END$/{//!p;}' $input > fragment2_EXTRA
sed -n '/^complex EXTRA$/,/^complex EXTRA END$/{//!p;}' $input > complex_EXTRA

submit="amspython \$HOSTPYFRAG/standalone/adf_new/PyFrag.py \\"

jobsubargue jobsub.txt                                      >> ./sub
echo $submit                                                >> ./sub
pyfragargue pyfrag.txt                                      >> ./sub



# Checks whether the old of new ADF input file is used, given by old_adfinputfile and adfinputfile, respectively 
Options=(fragment1_EXTRA fragment2_EXTRA complex_EXTRA adfinputfile old_adfinputfile)


for item in ${Options[*]}
do
  if [ -s $item ]; then
    subItem="--"${item}" $SCRIPTPATH/${item} \\"
    echo $subItem                                           >> ./sub
  else
    rm $item
  fi
done

rm jobsub.txt pyfrag.txt