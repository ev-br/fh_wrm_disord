#!/bin/bash
#SBATCH --job-name=m5.2L4
#SBATCH -n 4

#module load openblas
#mpirun ../../a.out.mpi

module load INTEL/parallel_studio_xe_2020_ce
mpirun ../../a.i.out.mpi

