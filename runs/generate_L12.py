from generate_runs import write_out_one, read_par_template

base_dct = {"L": 12, "amp": "0.1", "beta": "4.2",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 1, "stat": 0, "therm": "0",   # new or restart
            "step_p": "1d6", "step_w": "1d6",      # printout/checkpoint
            "time": "144.0",
       }

par_template = read_par_template("par_L12.template")

### beta = 4.2

sbatch_files = []
for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)
    

for sf in sbatch_files:
    print("sbatch " + sf)

### beta = 4.55
base_dct["beta"] = "4.55"

sbatch_files = []
for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)
    

for sf in sbatch_files:
    print("sbatch " + sf)

