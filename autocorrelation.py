# https://github.com/rauscher-lab/autocorrelation.git

import numpy as np
from scipy.signal import correlate


def one_hot_encode(num, max_num):
    vector = np.zeros(max_num, dtype=np.int32)
    if num > 0:
        vector[num-1] = True 
    return vector
    
def vector_autocorrelate(ts_array):
    n_timepoints = len(ts_array)
    n_dimensions = len(ts_array[0])
    acorr = [np.correlate(ts_array[:,i],ts_array[:,i],'full') for i in range(n_dimensions)]# correlate each component independently
    acorr = np.array(acorr)[:,n_timepoints-1:] # take the second half
    acorr = np.sum(acorr, axis = 0) # sum the correlations for each component
    acorr = acorr * 1./(n_timepoints - np.arange(n_timepoints)) # divide by the number of values actually measured and return
    return acorr

def vector_autocorrelate_signal(ts_array):
    n_time = ts_array.shape[0]
    n_dim = ts_array.shape[1]
    acorr = np.array([correlate(ts_array[:,i],ts_array[:,i],'full') for i in range(n_dim)])# correlate each component independently
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
    ac = vector_autocorrelate(oh_ts)
    return ac

def index_autocorrelate(mol_ts):
    
    unique_mols, indexes = np.unique(mol_ts, return_inverse=True)
    unique_indexes = np.unique(indexes)
        
    if np.all(unique_mols == 0) and len(unique_mols) == 1: # if no molecules found
        oh_ts = np.zeros((len(mol_ts), 1), dtype=np.int32)
        
    elif np.all(unique_indexes == 0) and len(unique_mols) == 1: # if the same one molecule found
        oh_ts = np.ones((len(mol_ts), 1), dtype=np.int32)
        
    elif np.all(unique_mols > 0) and  len(unique_mols) > 1: # if different with no zero
        oh_array = np.zeros((unique_indexes.shape[0] + 1, max((np.max(indexes + 1), 1))), dtype=np.int32)
        for un_index in unique_indexes + 1:
            oh_array[un_index] = one_hot_encode(un_index, max((np.max(indexes + 1), 1)))
        oh_ts = np.take(oh_array, indexes+1, axis=0)
        
    elif np.any(unique_mols == 0) and len(unique_mols) > 1: # if different ones (including 0)
        oh_array = np.zeros((unique_indexes.shape[0], max((np.max(indexes), 1))), dtype=np.int32)
        for un_index in unique_indexes:
            oh_array[un_index] = one_hot_encode(un_index, max((np.max(indexes), 1)))
        oh_ts = np.take(oh_array, indexes, axis=0) 
    ac = vector_autocorrelate_signal(oh_ts)
    return ac
