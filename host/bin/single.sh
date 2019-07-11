#!/bin/sh


JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

for prefile in adfinputfile sub
do
  if  [ -f $JOBDIR/$prefile ]; then
    rm -rf $JOBDIR/$prefile
  fi
done


bash "$HOSTPYFRAG"/standalone/adf_single/pyfragparce.sh "$1"

sbatch $JOBDIR/sub
