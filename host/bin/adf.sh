#!/bin/sh

JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

for prefile in adfinputfile sub complex_EXTRA fragment1_EXTRA fragment2_EXTRA
do
  if [ -f "$JOBDIR/$prefile" ]; then
    rm -rf "$JOBDIR/$prefile"
  fi
done

bash "$HOSTPYFRAG"/standalone/adf_new/pyfragparce.sh "$1"

sbatch $JOBDIR/sub || bash $JOBDIR/sub
