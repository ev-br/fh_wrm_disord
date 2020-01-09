#!/bin/bash
#SBATCH --job-name=m5.2L4
#SBATCH -s
#SBATCH -c 4
#SBATCH --partition=normal

#module load openblas

# The number of jobs MUST match the `SBATCH -c` # above

mpirun -n 4 ../../a.out.mpi

