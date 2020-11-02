#!/bin/bash
#SBATCH --job-name=L8b4.5r10a0.1
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L8b4.5r10a0.1
