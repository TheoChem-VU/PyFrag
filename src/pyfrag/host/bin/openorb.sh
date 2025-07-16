#!/bin/sh


JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

for prefile in sub
do
  if  [ -f $JOBDIR/$prefile ]; then
    rm -rf $JOBDIR/$prefile
  fi
done


bash "$HOSTPYFRAG"/src/pyfrag/host/standalone/adf_openorb/pyfragparce.sh "$1"

bash "$JOBDIR"/sub
