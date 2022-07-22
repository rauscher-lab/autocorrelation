# https://github.com/rauscher-lab/autocorrelation.git

import numpy as np
from scipy.signal import correlate
from scipy.optimize import curve_fit

def one_hot_encode(num, max_num):
    vector = np.zeros(max_num, dtype=np.bool)
    if num > 0:
        vector[num-1] = True 
    return vector
    
def vector_autocorrelate(ts_array, dtype=float):
    n_timepoints = len(ts_array)
    n_dimensions = len(ts_array[0])
    acorr = [np.correlate(ts_array[:,i].astype(dtype), ts_array[:,i].astype(dtype),'full') for i in range(n_dimensions)]# correlate each component independently
    acorr = np.array(acorr)[:,n_timepoints-1:] # take the second half
    acorr = np.sum(acorr, axis = 0) # sum the correlations for each component
    acorr = acorr * 1./(n_timepoints - np.arange(n_timepoints)) # divide by the number of values actually measured and return
    return acorr

def vector_autocorrelate_signal(ts_array, dtype=float):
    n_time = ts_array.shape[0]
    n_dim = ts_array.shape[1]
    acorr = np.array([correlate(ts_array[:,i].astype(dtype), ts_array[:,i].astype(dtype),'full') for i in range(n_dim)])# correlate each component independently
    acorr = acorr[:, n_time-1:]
    acorr = np.sum(acorr, axis = 0) # sum the correlations for each component
    acorr = acorr * 1./(n_time - np.arange(n_time)) # divide by the number of values actually measured and return
    return acorr

def index_autocorrelate_old(index_ts):
    it2id = {item: id_ for id_,item in enumerate(set(index_ts))}
    #id2it = {id_: item for id_,item in enumerate(set(index_ts))} 
    ts_corrected = np.array([it2id[item] for item in index_ts])
    oh_ts = [one_hot_encode(corrected_index, max(ts_corrected)) for corrected_index in ts_corrected]
    oh_ts = np.array(oh_ts)
    ac = vector_autocorrelate(oh_ts, dtype=np.int8)
    return ac

def index_autocorrelate(mol_ts):
    
    unique_mols, indexes = np.unique(mol_ts, return_inverse=True)
    unique_indexes = np.unique(indexes)
        
    if np.all(unique_mols == 0) and len(unique_mols) == 1: # if no molecules found
        oh_ts = np.zeros((len(mol_ts), 1), dtype=np.bool)
        
    elif np.all(unique_indexes == 0) and len(unique_mols) == 1: # if the same one molecule found
        oh_ts = np.ones((len(mol_ts), 1), dtype=np.bool)
        
    elif np.all(unique_mols > 0) and  len(unique_mols) > 1: # if different with no zero
        oh_array = np.zeros((unique_indexes.shape[0] + 1, max((np.max(indexes + 1), 1))), dtype=np.bool)
        for un_index in unique_indexes + 1:
            oh_array[un_index] = one_hot_encode(un_index, max((np.max(indexes + 1), 1)))
        oh_ts = np.take(oh_array, indexes+1, axis=0)
        
    elif np.any(unique_mols == 0) and len(unique_mols) > 1: # if different ones (including 0)
        oh_array = np.zeros((unique_indexes.shape[0], max((np.max(indexes), 1))), dtype=np.bool)
        for un_index in unique_indexes:
            oh_array[un_index] = one_hot_encode(un_index, max((np.max(indexes), 1)))
        oh_ts = np.take(oh_array, indexes, axis=0)
        
    ac = vector_autocorrelate_signal(oh_ts, dtype=float)
    return ac

def auto_window(taus, c):
    """
    This function estimates where to stop in summation of the estimator of AC time (Adapted from https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr)
    """
    m = np.arange(len(taus)) < c * taus
    if np.all(m):
        return -1
    else:
        return np.argmin(m)

def exp_f(t, tau):
    """The function to fit normalized AC with"""
    return np.exp(-t / tau)

def autocorrelation_time_sum(ac_norm, c=5):
    """
    Estimate autocorrelation time given AC function using a summation estimator as in (https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr)
    """
    taus = np.cumsum(ac_norm)
    windows = auto_window(taus, c)
    tau_est = taus[windows]
    return tau_est

def autocorrelation_time_fit(ac_norm):
    """
    Estimate AC time via fitting normalized AC to the exponential
    """
    popt, pcov = curve_fit(exp_f, np.arange(len(ac_norm)), ac_norm)
    tau_est = popt[0]
    return tau_est
