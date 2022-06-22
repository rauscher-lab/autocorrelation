# Autocorrelation

The `numpy.autocorrelate()` only works with one-dimensional time series. Sometimes we are interested in computing autocorrelation for time-series of multidimensional vectors. We can do that by considering each component independently and summing the autocorrelation over all components.  

This pipelines is designed to compute autocorrelation functions for time-dependent vectors. Make sure your vectors are normalized since vectors can be correlated in a sense of their direction. In addition to `vector_autocorrelate(ts_array)`, the code also contains the function `index_autocorrelate(index_ts)` which acts on the time-series of the indexes, for example, an atom or a molecule that visits a specific volume. 0 means nothing is present in the volume. The code can compute autocorrelation by one-hot-encoding all indexes and applying a simple vector correlation function.  The average lifetime of the same index (vector direction) can be estimated from the resulting autocorrelation by fitting it with exponential exp(-t/tau). tau is the average residence (life) time (of the *same* index in the volume).

## Usage:

Autocorrelation for time-series of normalized vectors:

```python
ts_array = np.array([
            [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.8, 0.6, 0.0], [0.0, 1.0, 0.0]
])
ac = vector_autocorrelate(ts_array)
print(ac) # [1, 0.8, 0.4, 0.0]
```

Autocorrelation for time-series of the indexes (0 means the molecules is absent from the volume):

```python
index_ts = [35, 35, 36, 0]
ac = index_autocorrelate(index_ts)
print(ac) # [3/4, 1/3, 0, 0]
```
The resulting `ac` shows the value of correlation of the signal with itself at different lagtime from `0` to `t`, if there are `t` time points in the array. For 3-dimensional vectors, `ac[0]==1` since the dot product of the normalized vector with itself equals to 1. However, for the case of indeces, `ac[0]` shows the occupancy (or, probability to find *any* index in the volume). 

## Added new features - estimating autocorrelation times
- Using an exponential fit of the normalized AC function (`autocorrelation_time_fit()`)
- Using an estimator based on AC series summation (`autocorrelation_time_sum()`) adapted from (https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr)

## Reading
https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Concepts_in_Biophysical_Chemistry_(Tokmakoff)/06%3A_Dynamics_and_Kinetics/22%3A_Biophysical_Reaction_Dynamics/22.05%3A_Time-Correlation_Functions
