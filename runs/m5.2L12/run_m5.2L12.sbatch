#!/bin/bash
#SBATCH --job-name=m5.2L12
#SBATCH -n 24
#SBATCH -t 5-00:00

#module load openblas
#mpirun ../../a.out.mpi

module load INTEL/parallel_studio_xe_2020_ce
mpirun ../../a.i.out.mpi

