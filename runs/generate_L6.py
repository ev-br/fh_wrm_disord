from generate_runs import write_out_one, read_par_template

base_dct = {"L": 6, "amp": "0.1", "beta": "3.5",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 1, "stat": 0, "therm": "0",   # new or restart
            "step_p": "5d7", "step_w": "5d7",      # printout/checkpoint
       }

par_template = read_par_template("par_L6.template")

sbatch_files = []
for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)
    
    
for sf in sbatch_files:
    print("sbatch " + sf)
