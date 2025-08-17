JOBSUB

#!/bin/bash
#SBATCH -J gausssian_pyfrag
#SBATCH -t 5:00
#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --partition=short
#SBATCH -o job-%j.stdout
#SBATCH -e job-%j.stderr

module load g09/A02

export NSCM=16

export WORK_TMPDIR=$TMPDIR
cd $WORK_TMPDIR

JOBSUB END

INPUT_SPECS
type = IRC
output file = Ethylene-forward.amv
fa1_name = complex
frag1 = C4H6
1.C
2.H
3.C
4.H
5.C
6.H
7.C
8.H
13.H
14.H
end frag1
frag2 = C2H4
9.C
10.C
11.H
12.H
15.H
16.H
end frag2

print bond 1 9 1.384
print strain frag1  1000
print strain frag2  2000

END INPUT_SPECS


"g09" <<eor


%nprocs=16
%mem=14000mb
#OPBE/6-31G*

Comments

0 1
END INPUT
