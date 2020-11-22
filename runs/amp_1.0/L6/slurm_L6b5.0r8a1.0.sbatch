#!/bin/bash
#SBATCH --job-name=L6b5.0r8a1.0
#SBATCH -n 1
#SBATCH -t 7-00:00

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L6b5.0r8a1.0
