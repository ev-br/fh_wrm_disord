"""
Generate a par_* and sbatch* scripts + the queueing script
"""
import os
import shutil
import stat


# @ HARISMA
ROOT = "~/fhdw/"

# @ duneyrr
ROOT = "~/sweethome/ferm/worm_disord/"

ROOT = os.path.expanduser(ROOT)

REPLICAS_STORE = os.path.join(ROOT, "replicas")
TARGET_PATH = os.path.join(ROOT, "runs/L4b3.5")


PAR_TEMPLATE = \
r"""__v2__    ! version tag
3         ! dimension
%(L)s         ! N(1)
%(L)s         ! N(2)
%(L)s         ! N(3)
%(cnf)s         ! 0 if new configuration, 1 if old one
%(stat)s         ! 0 if new statistics,    1 if old one
0         ! 0 if new rndm() seed
-5.2d0              ! \mu
%(amp)s       ! disorder amplitude
7.915d0  0.015d0     ! - U, initial for thermalization
%(beta)sd0               ! \beta
5.1d0      ! eta: GF vs Z weight
500       ! number of \tau points for GF tabulation
1.d2      ! tolerance level for determinant recalculation
0  0.1d0  ! x- and \tau- radii for cre/ann
1  0.1d0  ! x- and \tau- radii for leap_add/drop
20       ! nt
%(therm)s         ! number of sweeps for thermalization
2.d7     ! step for printing
6.d7      ! step for writing to disk
1.d0      ! step for measuring
6.0       ! time limit, hrs
------
0.0       ! add/drop
0.1       ! cre/ann
0.35      ! leap_add/drop
0.1     ! hop
------
%(seed1)s %(seed2)s
4836 2738
"""


SLURM_TEMPLATE = \
r"""#!/bin/bash
#SBATCH --job-name=%(suffix)s
#SBATCH -n 1

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_%(suffix)s
"""


def get_suffix(dct):
    return "L%(L)sb%(beta)sr%(replica)sa%(amp)s" % dct



base_dct = {"L": 4, "amp": "0.1", "beta": "3.5",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 0, "stat": 0, "therm": "1d2",   # new or restart
            "step_p": "2d7", "step_w": "6d7",      # printout/checkpoint
       }


#########################################
if __name__ == "__main__":

    sbatch_files = []
    
    for replica in range(1, 3):
        dct = base_dct.copy()
        dct["replica"] = replica
        dct["seed2"] += replica
        dct["suffix"] = get_suffix(dct)
        
        # write out the par file
        parfname = "par_%s" % dct["suffix"]
        parfname = os.path.join(TARGET_PATH, parfname)
        
        with open(parfname, "w") as parf:
            parf.write(PAR_TEMPLATE % dct)
        
        # write out the slurm file
        slurmfname = "slurm_%s.sbatch" % dct["suffix"]
        slurmfname = os.path.join(TARGET_PATH, slurmfname)
        with open(slurmfname, "w") as sf:
            sf.write(SLURM_TEMPLATE % dct)
        sbatch_files.append(slurmfname)
        
        # chmod u+x for SLURM
        st = os.stat(slurmfname)
        os.chmod(slurmfname, st.st_mode | stat.S_IEXEC)
        
        # copy the replica over
        r_path = "disord_L%sr%s.dat" % (dct["L"], replica)
        r_path = os.path.join(REPLICAS_STORE, r_path)
        shutil.copy(r_path,
                    os.path.join(TARGET_PATH, "disord_%s" % dct["suffix"]))

    # write out the run file to queue them all
    for sf in sbatch_files:
        print("sbatch " + sf)

