#!/bin/sh

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

cp $JOBDIR/sub $JOBDIR/sub_new

cat <<EOF >> "$JOBDIR/sub_new"
source "$HOSTPYFRAG"/src/pyfrag/venv_runner.sh && run_pyfrag_cluster "$HOSTPYFRAG"/src/pyfrag/cluster_runner.py --backend orca --job-dir "$JOBDIR" --input-file "$1"
EOF

sbatch $JOBDIR/sub_new

rm $JOBDIR/sub_new
