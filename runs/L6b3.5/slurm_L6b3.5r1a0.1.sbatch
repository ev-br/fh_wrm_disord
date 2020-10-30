#!/bin/bash
#SBATCH --job-name=L6b3.5r1a0.1
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L6b3.5r1a0.1
