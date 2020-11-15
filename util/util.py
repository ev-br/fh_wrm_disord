import io

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


__all__ = ["Bunch", "fit_Tc"]

########## Fit Tc : Goulko, Wingate, arXiv:1008.3348

class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def fitfunc_Tc(LT, f0, f1, c, Tc):
    """Eq (10) of arXiv:1008.3348."""
    w = 0.8
    nu = 0.67
    nu1 = 1. / nu
    
    L, T = LT
    res = f0 + f1*(T - Tc)*L**nu1
    res *= 1. + c*L**(-w)
    return res


def fit_Tc(df):
    popt, pcov = curve_fit(fitfunc_Tc, (df["L"], df["T"]), df["g_im"],
                       sigma=df["err_g_im"], absolute_sigma=True)
    
    Tc = popt[-1]
    err_Tc = np.sqrt(pcov[-1, -1])
    
    # beta = 1 / T   ==>  err_beta = err_T / T**2
    beta_c = 1. / Tc
    err_beta_c = err_Tc / Tc**2
    
    return Bunch(Tc=Tc, err_Tc=err_Tc,
                 beta_c=beta_c, err_beta_c=err_beta_c,
                 popt=popt, pcov=pcov)


def replica_average(mean, err, num_samples=100000, seed=1234):
    """Compute the replica average over a `mean` array w/synthetic datasets.
        
    Given two arrays of values and errorbars (one for each replica),
    generate `num_samples` synthetic realizations for each replica. For each realization,
    compute the average value.
    
    Parameters
    ----------
    mean : array_like, shape (n_meas,)
        Mean values, one for each replica.
    err : array_like, shape (n_meas,)
        Errorbars, one for each replica.
    num_samples : int
        Generate sythetic datasets of this many samples
        for each element of `mean` and arr.
    seed : int
        np.random seed for random numbers
        
    Returns
    -------
    ndarray, shape (n_samples,)
        Synthetic dataset: average over the replicas
        
    """
    mean, err = map(np.asarray, (mean, err))
    num_meas = mean.size
    synth = np.empty((num_meas, num_samples), dtype=float)
    
    rndm = np.random.RandomState(seed)
    for j in range(num_meas):
        synth[j, :] = rndm.normal(mean[j], err[j], size=num_samples)
        
    return synth.sum(axis=0) / num_meas


############## Process MC data for the worm/disord code

def read_res(fnames):
    """Read the `res_SUFFIX` files for a collection of replicas."""
    
    firstline = "replica_id density err_density conv_density g_im err_g_im conv_g_im Z(mln)"
    lines = [firstline]

    for fname in fnames:
        with open(fname, 'r') as f:
            line = f.read()
        lines.append(" ".join(line.split()))

    buf = io.StringIO("\n".join(lines))
    df = pd.read_csv(buf, sep="\s+")
    return df




