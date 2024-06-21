#!/bin/sh


JOBDIR="$( cd "$(dirname "$1")" ; pwd -P )"

for prefile in adfinputfile sub complex_EXTRA fragment1_EXTRA fragment2_EXTRA
do
  if  [ -f $JOBDIR/$prefile ]; then
    rm -rf $JOBDIR/$prefile
  fi
done


bash "$HOSTPYFRAG"/standalone/adf_newopen/pyfragparce.sh "$1"

# The user may run this file on a cluster (sbatch) or a local machine (bash).
# However, the user must have either sbatch or bash installed (which may not be the case on Windows computers).
# if command -v sbatch >/dev/null 2>&1; then
sbatch "$JOBDIR/sub"
# elif command -v bash >/dev/null 2>&1; then
#     bash "$JOBDIR/sub"
# else
#     echo "Warning: Please activate a shell that supports sbatch or bash."
# fi