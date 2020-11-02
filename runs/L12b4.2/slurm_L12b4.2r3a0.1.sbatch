#!/bin/bash
#SBATCH --job-name=L12b4.2r3a0.1
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_L12b4.2r3a0.1
