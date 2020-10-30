from generate_runs import write_out_one

base_dct = {"L": 6, "amp": "0.1", "beta": "3.5",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 0, "stat": 0, "therm": "1e2",   # new or restart
            "step_p": "5d7", "step_w": "5d7",      # printout/checkpoint
       }

sbatch_files = []
for replica in range(1, 11):
    sf = write_out_one(base_dct, replica)
    sbatch_files.append(sf)
    
    
for sf in sbatch_files:
    print("sbatch " + sf)
