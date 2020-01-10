#!/bin/bash
#SBATCH --job-name=m5.2L4
#SBATCH -n 4

module load openblas
mpirun ../../a.out.mpi

