from generate_runs import write_out_one, read_par_template

base_dct = {"L": 8, "amp": "0.1", "beta": "4.0",   # physical
            "seed1": 4836, "seed2": 2738,
            "replica" : "1",
            "cnf": 1, "stat": 0, "therm": "0",   # new or restart
            "step_p": "2d6", "step_w": "2d6",      # printout/checkpoint
       }

par_template = read_par_template("par_L8.template")

sbatch_files = []
for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)

    
### beta = 4.25
base_dct["beta"] = "4.25"

for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)


### beta = 4.5
base_dct["beta"] = "4.5"

for replica in range(1, 11):
    sf = write_out_one(base_dct, replica, par_template)
    sbatch_files.append(sf)


##############################
for sf in sbatch_files:
    print("sbatch " + sf)

