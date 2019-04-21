#!/bin/sh

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

cp $JOBDIR/sub $JOBDIR/sub_new

cat <<EOF >> "$JOBDIR/sub_new"
python3 "$HOSTPYFRAG"/standalone/gaussian/pyfrag.py "$JOBDIR"/$1
EOF

sbatch $JOBDIR/sub_new

rm $JOBDIR/sub_new
