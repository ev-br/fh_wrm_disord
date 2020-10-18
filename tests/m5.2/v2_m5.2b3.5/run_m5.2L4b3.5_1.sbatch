#!/bin/bash
#SBATCH --job-name=m5.2L4b3.5_1
#SBATCH -n 1

#module load openblas
#mpirun ../../a.out.mpi

module load INTEL/parallel_studio_xe_2020_ce
#mpirun ../../a.i.out.mpi

~/fhdw/a.out par_m5.2L4b3.5_1

