Special Pyfrag Calculation
==========================

Besides the above simple calculations, it is more complicated to perform an open shell Activation Strain Analysis (ASA) using PyFrag 2019 for the technical reasons. For more information please check the example_ consisting of an analysis of the C-C single bond between two CP radicals in the four-atomic molecule PCCP.

Open Shell ASA
--------------

The basic PyFrag 2019 input for the Activation Strain Analysis (ASA) using ADF to perform an open shell Activation Strain Analysis is as follow: ::

   JOBSUB

   #!/bin/bash
   #SBATCH -J NNC
   #SBATCH -N 1
   #SBATCH -t 24:00:00
   #SBATCH --ntasks-per-node=24
   #SBATCH --partition=normal
   #SBATCH --output=%job.stdout
   #SBATCH --error=%job.stdout
   export NSCM=24

   JOBSUB END


   PyFrag

   ircpath /Users/xiaobo/Desktop/test/molecule.xyz
   fragment  1 3 4 5
   fragment  2 6 7 8
   strain    0
   strain    0
   bondlength 1 2  1.52

   PyFrag END

   fragment1 EXTRA

   SYMMETRY C(3V)
   CHARGE    0 0

   OCCUPATIONS
   E1 4
   A1 5
   END

   fragment1 EXTRA END

   fragment2 EXTRA

   SYMMETRY C(3V)
   CHARGE    0 0

   OCCUPATIONS
   E1 4
   A1 5
   END

   fragment2 EXTRA END

   complex EXTRA

   FRAGOCCUPATIONS

   frag1
   E1 2//2
   A1 3//2
   SUBEND

   frag2
   E1 2//2
   A1 2//3
   SUBEND

   END

   complex EXTRA END

   fragment1 open EXTRA
   charge 0 1
   unrestricted
   fragment1 open EXTRA END

   fragment2 open EXTRA
   charge 0 1
   unrestricted
   fragment2 open EXTRA END

   complex open EXTRA
   charge 0 0
   complex open EXTRA END


   ADF

   basis
   type DZ
   core None
   end

   xc
   gga OPBE
   end


   scf
   iterations 99
   converge 0.0001
   mixing 0.20
   end

   numericalquality good


   ADF END


To submit a job, create a directory and generate a input file and run the following command to submit a job:

``pyfrag -x open job.in``


In order to perform a successful open shell fragment analysis, additional information should be provided in the following input blocks: ::

   fragment1 EXTRA

   SYMMETRY C(3V)
   CHARGE    0 0

   OCCUPATIONS
   E1 4
   A1 5
   END

   fragment1 EXTRA END

   fragment2 EXTRA

   SYMMETRY C(3V)
   CHARGE    0 0

   OCCUPATIONS
   E1 4
   A1 5
   END

   fragment2 EXTRA END

   complex EXTRA

   FRAGOCCUPATIONS

   frag1
   E1 2//2
   A1 3//2
   SUBEND

   frag2
   E1 2//2
   A1 2//3
   SUBEND

   END
   complex EXTRA END


The fragment calculations used to provide the TAPE21 for the overall complex calculation
must be done, for technical reasons, in the restricted mode. The proper spins are then
specified in the calculation of the overall molecule using the FragOccupations key.
Noted a proper decomposition of an electron-pair bond energy requires specifying
opposite spins for the unpaired electrons of the respective radical fragments,
which can be done with the input key FragOccupations. For the convenience of the analysis,
it is suggested to specify the electronic configuration according to the symmtry of the molecule.

Please note that if one neglects explicitly specifying opposite spins for the
unpaired electrons of the fragments, each of them is treated as being half an
alpha and half a beta electron and consequently, they enter into a spurious
Pauli repulsive interaction. This results, among others, into the Pauli
repulsion term being too repulsive and the orbital interaction term being too
much stabilizing.

Note that this implies a slight approximation because the bond energy
computed in this way refers to the energy difference between complex
and two fragment radicals that are described by orbitals from a spin-restricted SCF
calculation, which have been given an unrestricted occupation. In other words, the
set of alpha- and beta-spin orbitals are identical and the effect of spin polarization is missing.
In practice, this leads to minor energy differences with respect to the correct
bond energy, that is, the energy difference between complex and two
fragment radicals treated in the unrestricted mode, i.e., for which the set of
alpha- and beta-spin orbitals are allowed to relax toward different solutions
in the SCF procedure.

This correction term can be computed directly by carrying out an unrestricted computation of the
fragment radical using the following block: ::

   fragment1 open EXTRA
   charge 0 1
   unrestricted
   fragment1 open EXTRA END

   fragment2 open EXTRA
   charge 0 1
   unrestricted
   fragment2 open EXTRA END

   complex open EXTRA
   charge 0 0
   complex open EXTRA END


After the calculation, all results will be summarized in two text files. One file with the name started with pyfrag1
include all terms obtained from the above open shell ASA.

The second file with the name started with pyfrag2 include the correction energy terms from the correction procedure later.


New Open Shell ASA (Since ADF 2019)
-----------------------------------

Since ADF 2019, new method to do open-shell fragment analysis has been included. For details, please refer to the ADF_ website.

Based on this method, a new module to do the activation strain analysis has been developed by Xiaobo Sun and Eva Blokker. The specification for the print options is similar with the previous one, except to has to specify the spin state of orbital, such as 1_A, which means spin-A orbital 1. All the following options are acceptable: ::

   overlap frag1 HOMO frag2 HOMO
   overlap A1 frag1 3_B S frag2 1_B

   orbitalenergy frag1 HOMO-2
   orbitalenergy frag1 HOMO-1
   orbitalenergy frag1 LUMO
   orbitalenergy frag2 LUMO
   orbitalenergy A1 frag1 3_A

   population frag1 HOMO
   population frag2 HOMO
   population frag2 LUMO
   population A1 frag1 3_A
   population S frag2 1_B


The basic PyFrag 2019 input for the Activation Strain Analysis (ASA) using ADF 2019 to perform a new open shell Activation Strain Analysis is as follow: ::


   JOBSUB
   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH --ntasks-per-node=16
   #SBATCH --partition=tc
   #SBATCH --time=24:00:00
   #SBATCH --job-name=methane
   #SBATCH --output=methane.out
   #SBATCH --error=methane.err
   module load adf/2019.301
   JOBSUB END

   ADF
   XC
     GGA BLYP
     DISPERSION Grimme3 BJDAMP
   END

   NumericalQuality Excellent

   BASIS
     TYPE TZ2P
     CORE None
   END

   SCF
     ITERATIONS 300
   END

   SYMMETRY AUTO
   CHARGE 0
   ADF END

   PyFrag
   ircpath /home/x2sun/methane.amv
   fragment  1 2 3 4
   fragment  5
   strain    0
   strain    0
   bondlength 1 5
   overlap A1 frag1 3_A S frag2 1_A
   overlap A1 frag1 3_B S frag2 1_B
   overlap A1 frag1 2_B S frag2 1_B
   population frag1 HOMO
   population frag2 HOMO
   PyFrag END


   fragment1 EXTRA
   SYMMETRY C(3V)
   CHARGE 0 1
   unrestricted

   IrrepOccupations
   E1 2//2
   A1 3//2
   END
   fragment1 EXTRA END


   fragment2 EXTRA
   SYMMETRY AUTO
   CHARGE 0 -1
   Unrestricted

   IrrepOccupations
   S 0//1
   END
   fragment2 EXTRA END


   complex EXTRA
   UnrestrictedFragments
   unrestricted
   complex EXTRA END

The molecule is methane: ::

   C      -0.88533700      -1.60854000       0.00000000
   H      -0.50220300      -2.11092900       0.89352900
   H      -0.50220300      -2.11092900      -0.89352900
   H      -1.97897100      -1.64799500       0.00000000
   H      -0.55799500      -0.56431300       0.00000000


To submit a job, create a directory and generate a input file and run the following command to submit a job:

``pyfrag -x newopen job.in``



Open Shell ASA Orbital Energy
-----------------------------

Because the above open shell Activation Strain Analysis will not give the correct orbital energy of fragment,
thus, in order to extract the correct orbital energy, the following small calculation can be performed: ::


   PyFrag

   ircpath /Users/xiaobo/Desktop/test/plams.0001
   fragment  frag1open
   orbitalenergy  HOMO
   orbitalenergy  HOMO-1
   orbitalenergy  LUMO
   orbitalenergy  LUMO+1
   orbitalenergy  AA 5

   PyFrag END

The ircpath refer to the plams directory that contains all the open shell calculation results.
Besides, the fragment term specifies from which (fragment1open or fragment1open) orbital energy will be extracted.
Noted only one fragment informatin can be readed for one calculation.


To submit a job, create a directory and generate a input file and run the following command to run a job:

``pyfrag -x openorb job.in``


Single Points
-------------

The basic PyFrag 2019 input for the Activation Strain Analysis (ASA) using ADF to do single point calculation for a series of coordinates is as follows: ::

   JOBSUB

   #!/bin/bash
   #SBATCH -J NNC
   #SBATCH -N 1
   #SBATCH -t 1:00:00
   #SBATCH --ntasks-per-node=24
   #SBATCH --partition=short
   #SBATCH --output=%job.stdout
   #SBATCH --error=%job.stdout
   export NSCM=24

   JOBSUB END

   PyFrag

   ircpath /Users/xiaobo/Desktop/test1/molecule.xyz

   VDD 1 2 3
   angle 1 2 3 90
   bondlength 1 2 5

   PyFrag END


   ADF

   basis
   type DZ
   core None
   end

   xc
   gga OPBE
   end

   scf
   iterations 99
   converge 0.0001
   mixing 0.20
   end

   numericalquality good

   ADF END


Note that the fragment definations are not needed. This functionality provide an easy way to do a simple single point calculation for a series of different molecular coordinates and get the computational results like VDD charges, total energy, bond length and angles. Use the following command to run this calculation:

``pyfrag -x single job.in``


.. _example: https://www.scm.com/doc/ADF/Examples/PCCP_Unr_BondEnergy.html?highlight=open+shell+fragment
.. _ADF: https://www.scm.com/doc/ADF/Examples/EDA_Unr_CH3I.html#example-eda-unr-ch3i
