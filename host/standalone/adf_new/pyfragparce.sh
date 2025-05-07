### General information ###
# This script is used create an parser input file for the PyFrag.py script
# The pyfrag input file (*not* the parser input file called "sub") is split into three main sections
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
      # Skip lines that are fully commented (but allow lines starting with #SBATCH or #!/) (e.g., shebang)
      if [[ "$line" =~ ^# && ! "$line" =~ ^#SBATCH && ! "$line" =~ ^#!/ ]]; then
          continue
      fi

      # Remove inline comments (e.g., #, ;, or /) but preserve the line if it starts with #SBATCH or #!/  (e.g., shebang)
      if [[ ! "$line" =~ ^#SBATCH && ! "$line" =~ ^#!/ ]]; then
          line=$(echo "$line" | sed -E 's/[#;/].*$//' | sed '/^\s*$/d')
      fi

      # Process non-empty lines
      if [ "$line" != "" ]; then
          name="$line"
          echo "$name"
      fi
  done < "$jobsub"
}


function pyfragargue {
    # This function translates the PyFrag sections into a string of arguments (format: "--parser_key") that ends up in the "sub" file
    # Lines starting with # or ; are ignored, and inline comments are removed only if separated by a space
    pyfrag="$1"
    while read -r line
    do
        # Remove inline comments only if separated by a space
        line=$(echo "$line" | sed -E 's/\s[#;/].*$//' | sed '/^\s*$/d')

        # Skip empty lines
        if [ "$line" == "" ]; then
            continue
        fi

        name="$line"
        echo "--$name \\"

    done < "$pyfrag"
}

input=$*
SCRIPTPATH="$( cd "$(dirname "$1")" ; pwd -P )"

# Extract the relevant sections from the input file, e.g. JOBSUB, AMS, ADF, PyFrag, fragment1 EXTRA, fragment2 EXTRA, complex EXTRA
# Each section has a header (one line) and a footer (specified with END block)
# The sections are read in a case-insensitive manner through the sed command with the I flag
sed -n '/^JOBSUB$/I,/^JOBSUB END$/I{//!p;}' $input > jobsub.txt
sed -n '/^AMS$/I,/^AMS END$/I{//!p;}' $input > adfinputfile
sed -n '/^ADF$/I,/^ADF END$/I{//!p;}' $input > old_adfinputfile
sed -n '/^PyFrag$/I,/^PyFrag END$/I{//!p;}' $input > pyfrag.txt
sed -n '/^fragment1 EXTRA$/I,/^fragment1 EXTRA END$/I{//!p;}' $input > fragment1_EXTRA
sed -n '/^fragment2 EXTRA$/I,/^fragment2 EXTRA END$/I{//!p;}' $input > fragment2_EXTRA
sed -n '/^complex EXTRA$/I,/^complex EXTRA END$/I{//!p;}' $input > complex_EXTRA

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