import numpy as np


if __name__ == "__main__":

    L = 12
    seed = 1

    for seed in range(1, 11):
        rndm = np.random.RandomState(seed)
        replica = rndm.uniform(-1, 1, size=L*L*L)

        with open("disord_L%sr%s.dat" % (L, seed), "w") as f:
            f.write(str(seed)+"\n")
            f.write(str(L)+"\n")
            for value in replica:
                f.write(str(value)+"  ")
