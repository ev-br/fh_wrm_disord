#!/bin/bash
#SBATCH --job-name=m5.2L8
#SBATCH -n 32
#SBATCH -t 3-03:00

#module load openblas
#mpirun ../../a.out.mpi

module load INTEL/parallel_studio_xe_2020_ce
mpirun ../../a.i.out.mpi

