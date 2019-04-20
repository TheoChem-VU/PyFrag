Basic Usage Tutorial
====================

Usage
-----

The user can type pyfrag -h to see all the commands that can be used in this program, which will show: ::

   Usage: pyfrag [-h] [-s] [-x command]  [...]
   -h          : print this information
   -s          : run job quietly
   -x          : start the executable named command
               : command include restart, which restart job
               : restart, which restart a job after it is stoped
               : summary, which summarize all job result after jobs finished
               : default command is pyfrag itself
   The example command is like as follow, in which job.in is job input
   pyfrag job.in
   or
   pyfrag -x restart job.in
   or
   pyfrag -s -x summary job.in

To submit a job, create a directory and generate a input file and run the following command to submit a job. Note for each job, a new directory and a new job name should be given. Note: The user should avoid running more than one job in a single directory.

``pyfrag job.in``

To obtain the latest information about your job, the user can run:

``pyfrag -x summary job.in``

If a change in the input file is required, make the change and the resubmit the job using:

``pyfrag -x restart job.in``


Sample Input Example
--------------------

A simple job input is provided below. The input script can be roughly divided into four section: the required submit information for a job scheduling system (Slurm in this example), ADF parameters, PyFrag 2019 parameters, and geometry parameters. ::

   ''''
   JOBSUB section is for the information passed to the remote host machine
   where the heavy computational job is done! It is written in the fashion of Slurm.
   ''''
   JOBSUB

   #!/bin/bash
   #SBATCH -J frag_1
   #SBATCH -N 1
   #SBATCH -t 50:00
   #SBATCH --ntasks-per-node=24
   #SBATCH --partition=short
   #SBATCH --output=%job.stdout
   #SBATCH --error=%job.stdout
   export NSCM=24

   JOBSUB END

   ''''
   Provide the parameters for a DFT calculation using ADF software.
   ''''
   ADF

   basis
   type TZ2P
   core Small
   end

   xc
   gga OPBE
   end

   relativistic SCALAR ZORA

   scf
   iterations 299
   converge 0.00001
   mixing 0.20
   end

   numericalquality verygood

   charge 0 0
   symmetry auto


   ADF END

   ''''
   Provide the parameters for an activation strain analysis.
   Noted a bondlength calculation is needed to provilde x axis value for ASA.
   ''''

   PyFrag

   fragment  2
   fragment  1 3 4 5 6
   strain    0
   strain   -554.09
   bondlength 1 6  1.09

   PyFrag END


   ''''
   Guessed geometry coordinate for reactent1, reactent2, reactent complex,
   transition state and product.
   ''''

   Geometrycoor

   R1: Fe-II(CO)4 + CH4
   Pd       0.00000000       0.00000000       0.32205546

   R2: CH4
   C       0.00000000       0.00000000      -1.93543634
   H      -0.96181082       0.00000000      -1.33610429
   H       0.00000000      -0.90063254      -2.55201285
   H       0.00000000       0.90063254      -2.55201285
   H       0.96181082       0.00000000      -1.33610429

   RC: Fe-II(CO)4 + CH4
   C       0.00000000       0.00000000      -1.93543615
   Pd       0.00000000       0.00000000       0.322055
   H      -0.96181082       0.00000000      -1.33610429
   H       0.00000000      -0.90063254      -2.55201285
   H       0.00000000       0.90063254      -2.55201285
   H       0.96181082       0.00000000      -1.33610429

   TS: Fe-II(CO)4 + CH4
   C      -1.74196777      -2.22087997       0.00000000
   Pd     -2.13750904      -0.23784341       0.00000000
   H      -2.80956968      -2.49954731       0.00000000
   H      -1.26528821      -2.62993236       0.8956767
   H      -1.26528821      -2.62993236      -0.895676
   H      -0.75509932      -0.88569836       0.00000000

   P: Fe-II(CO)4 + CH4
   C      -2.10134690      -2.41901732       0.1862099
   Pd      -2.73145901      -0.57025833       0.419766
   H      -3.88639130      -1.04648079      -0.43099501
   H      -2.78392696      -3.12497645       0.66994616
   H      -1.97386865      -2.66955518      -0.87144525
   H      -1.12556673      -2.41201402       0.698583

   Geometrycoor END

The user  might want to specify an additional input for the different sections of the overall workflow. To specify additional information for say, fragment1 and fragment2 see the syntax shown below. Additional complex insert statements for the fragment analysis calculation can be added. Similarly, the R1 EXTRA, R2 EXTRA, RC EXTRA, TS EXTRA, P EXTRA, IR EXTRA insert statements for R1, R2, RC, TS, P, IRC calculation. ::

   fragment1 EXTRA
   charge 1
   fragment1 EXTRA END

   fragment2 EXTRA
   charge -1
   fragment2 EXTRA END

   complex EXTRA
   charge 2
   complex EXTRA END

   R1 EXTRA
   charge 0
   R1 EXTRA END

   R2 EXTRA
   charge 0
   R2 EXTRA END

   RC EXTRA
   charge 0
   RC EXTRA END

   TS EXTRA
   charge 0
   tsrc
   Bond 1 2 -1
   end
   TS EXTRA END

   P EXTRA
   charge 0
   P EXTRA END

   IR EXTRA
   Geometry
    IRC Backward POINTS=20 STEP=1
   ITERATIONS 300
   CONVERGE 0.000001
   End
   IR EXTRA END


Result example
--------------
After the job has been submitted, a website as shown in the figure below will be launched. The website summarizes all relevant information, including: a) the convergence criteria, b) the latest structure from the optimization in the form of movie, c) the latest energy and coordinates, and d) the activation strain analysis (once the complete workflow has finished). The user can decide if the optimization process is correct or incorrect, and if necessary, can stop the  job. If the input file is then modified or updated, the job will be resubmitted and the overall workflow will resume from where it left off.

.. image:: jobresult.png
   :alt: result
