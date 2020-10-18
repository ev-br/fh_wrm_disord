#!/bin/bash
#SBATCH --job-name=m5.2L4b3.5_2
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_m5.2L4b3.5_2

