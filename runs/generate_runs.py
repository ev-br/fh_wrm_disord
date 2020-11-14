"""
Generate a par_* and sbatch* scripts + the queueing script
"""
import os
import shutil
import stat


# @ HARISMA
ROOT = "~/fhdw/"

# @ duneyrr
#ROOT = "~/sweethome/ferm/worm_disord/"

ROOT = os.path.expanduser(ROOT)


SLURM_TEMPLATE = \
r"""#!/bin/bash
#SBATCH --job-name=%(suffix)s
#SBATCH -n 1
#SBATCH -t 7-00:00

module load INTEL/parallel_studio_xe_2020_ce

~/fhdw/a.out par_%(suffix)s
"""


def get_suffix(dct):
    return "L%(L)sb%(beta)sr%(replica)sa%(amp)s" % dct



base_dct = {"L": 4, "amp": "0.1", "beta": "3.5",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 0, "stat": 0, "therm": "1e2",   # new or restart
            "step_p": "2d7", "step_w": "6d7",      # printout/checkpoint
       }



def write_out_one(base_dct, replica, par_template):
    """Write out a single par/slurm file, return the sbatch line."""

    dct = base_dct.copy()
    dct["replica"] = replica
    dct["seed2"] += replica
    dct["suffix"] = get_suffix(dct)

    replicas_store = os.path.join(ROOT, "replicas")
    target_path = os.path.join(ROOT, "runs/L%(L)sb%(beta)s" % dct)

    # write out the par file
    parfname = "par_%s" % dct["suffix"]
    parfname = os.path.join(target_path, parfname)
    
    with open(parfname, "w") as parf:
        parf.write(par_template % dct)

    # sanity check:
    # FIXME
    if float(dct["stat"]) != 0 and float(dct["therm"]) != 0:
        print("!!! ", parfname, "stat & therm")

    # write out the slurm file
    slurmfname = "slurm_%s.sbatch" % dct["suffix"]
    slurmfname = os.path.join(target_path, slurmfname)
    with open(slurmfname, "w") as sf:
        sf.write(SLURM_TEMPLATE % dct)
    
    # chmod u+x for SLURM
    st = os.stat(slurmfname)
    os.chmod(slurmfname, st.st_mode | stat.S_IEXEC)
    
    # copy the replica over
    r_path = "disord_L%sr%s.dat" % (dct["L"], replica)
    r_path = os.path.join(replicas_store, r_path)
    shutil.copy(r_path,
                os.path.join(target_path, "disord_%s.dat" % dct["suffix"]))
    # FIXME: create a replica if not exists

    return slurmfname


def read_par_template(fname):
    with open(fname, 'r') as f:
        template = f.read()
    return template


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
                    os.path.join(TARGET_PATH, "disord_%s.dat" % dct["suffix"]))

    # write out the run file to queue them all
    for sf in sbatch_files:
        print("sbatch " + sf)

