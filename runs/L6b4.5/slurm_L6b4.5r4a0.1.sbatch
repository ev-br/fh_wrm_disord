#!/bin/bash
#SBATCH --job-name=L6b4.5r4a0.1
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L6b4.5r4a0.1
