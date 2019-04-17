#!/bin/sh

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

cp $JOBDIR/sub $JOBDIR/sub_new

cat <<EOF >> "$JOBDIR/sub"
python3 "$HOSTPYFRAG"/standalone/turbomole/pyfrag.py $1
EOF

sbatch $JOBDIR/sub_new

rm $JOBDIR/sub_new

