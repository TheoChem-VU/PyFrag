#!/bin/sh


JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

bash "$HOSTPYFRAG"/standalone/adf_new/pyfragparce.sh "$1"

sbatch $JOBDIR/sub
