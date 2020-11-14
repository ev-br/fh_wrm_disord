#!/bin/bash
#SBATCH --job-name=L6b4.5r5a1.0
#SBATCH -n 1
#SBATCH -t 3-00:00

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L6b4.5r5a1.0
