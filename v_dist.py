#!/usr/bin/env python

"""Calculate v baselines distribution.

:Authors: Shifan Zuo
:Date: 2015-01-08
:email: sfzuo@bao.ac.cn
"""

import argparse
import numpy as np
import h5py
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt


def unique(ar, return_index=False, return_inverse=False, return_counts=False):
    """
    Find the unique elements of an array.

    Returns the sorted unique elements of an array. There are two optional
    outputs in addition to the unique elements: the indices of the input array
    that give the unique values, and the indices of the unique array that
    reconstruct the input array.

    Parameters
    ----------
    ar : array_like
        Input array. This will be flattened if it is not already 1-D.
    return_index : bool, optional
        If True, also return the indices of `ar` that result in the unique
        array.
    return_inverse : bool, optional
        If True, also return the indices of the unique array that can be used
        to reconstruct `ar`.
    return_counts : bool, optional
        .. versionadded:: 1.9.0
        If True, also return the number of times each unique value comes up
        in `ar`.

    Returns
    -------
    unique : ndarray
        The sorted unique values.
    unique_indices : ndarray, optional
        The indices of the first occurrences of the unique values in the
        (flattened) original array. Only provided if `return_index` is True.
    unique_inverse : ndarray, optional
        The indices to reconstruct the (flattened) original array from the
        unique array. Only provided if `return_inverse` is True.
    unique_counts : ndarray, optional
        .. versionadded:: 1.9.0
        The number of times each of the unique values comes up in the
        original array. Only provided if `return_counts` is True.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> np.unique([1, 1, 2, 2, 3, 3])
    array([1, 2, 3])
    >>> a = np.array([[1, 1], [2, 3]])
    >>> np.unique(a)
    array([1, 2, 3])

    Return the indices of the original array that give the unique values:

    >>> a = np.array(['a', 'b', 'b', 'c', 'a'])
    >>> u, indices = np.unique(a, return_index=True)
    >>> u
    array(['a', 'b', 'c'],
           dtype='|S1')
    >>> indices
    array([0, 1, 3])
    >>> a[indices]
    array(['a', 'b', 'c'],
           dtype='|S1')

    Reconstruct the input array from the unique values:

    >>> a = np.array([1, 2, 6, 4, 2, 3, 2])
    >>> u, indices = np.unique(a, return_inverse=True)
    >>> u
    array([1, 2, 3, 4, 6])
    >>> indices
    array([0, 1, 4, 3, 1, 2, 1])
    >>> u[indices]
    array([1, 2, 6, 4, 2, 3, 2])

    """
    ar = np.asanyarray(ar).flatten()

    optional_indices = return_index or return_inverse
    optional_returns = optional_indices or return_counts

    if ar.size == 0:
        if not optional_returns:
            ret = ar
        else:
            ret = (ar,)
            if return_index:
                ret += (np.empty(0, np.bool),)
            if return_inverse:
                ret += (np.empty(0, np.bool),)
            if return_counts:
                ret += (np.empty(0, np.intp),)
        return ret

    if optional_indices:
        perm = ar.argsort(kind='mergesort' if return_index else 'quicksort')
        aux = ar[perm]
    else:
        ar.sort()
        aux = ar
    flag = np.concatenate(([True], aux[1:] != aux[:-1]))

    if not optional_returns:
        ret = aux[flag]
    else:
        ret = (aux[flag],)
        if return_index:
            ret += (perm[flag],)
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            iperm = perm.argsort()
            ret += (np.take(iflag, iperm),)
        if return_counts:
            idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
            ret += (np.diff(idx),)
    return ret



conf = ['3x32_0.2', '3x32_0.4', '3x32_1.0', '31_32_33_0.4', '3x32_0.4_alt', '3x32_81011']
Nconf = len(conf)
cyl1_fp = [] # feed v position for feed 1, unit: m
cyl2_fp = []
cyl3_fp = []


# Read in arguments
parser = argparse.ArgumentParser(description="Calculate v baselines distribution.")
parser.add_argument('-c', '--case', type=int, choices=range(1, Nconf+1), default=1, help='Which array configuration.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-f', '--freq', type=float, default=750.0, help='Observing frequency.')
parser.add_argument('-o', '--outfile', help='Output name file name.')
args = parser.parse_args()


c = 3.0e8 # light speed, m/s
wavelen = c / (args.freq * 1.0e6) # m, wavelength at 750.0 MHz

# Case 1: 3x32 D = 0.2
if args.case == 1:
    Nfeeds = 32
    D = 0.2 # m
    Dl = D / wavelen
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * Dl)
        cyl2_fp.append(fd * Dl)
        cyl3_fp.append(fd * Dl)
# Case 2: 3x32 D = 0.4
elif args.case == 2:
    Nfeeds = 32
    D = 0.4 # m
    Dl = D / wavelen
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * Dl)
        cyl2_fp.append(fd * Dl)
        cyl3_fp.append(fd * Dl)
# Case 3: 3x32 D = 1.0
elif args.case == 3:
    Nfeeds = 32
    D = 1.0 # m
    Dl = D / wavelen
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * Dl)
        cyl2_fp.append(fd * Dl)
        cyl3_fp.append(fd * Dl)
# Case 4: 31+32+33 D = 0.4
elif args.case == 4:
    Nfeeds1 = 31
    Nfeeds2 = 32
    Nfeeds3 = 33
    D1 = (Nfeeds2 - 1) * 0.4 / (Nfeeds1 - 1) # m
    D2 = 0.4 # m
    D3 = (Nfeeds2 - 1) * 0.4 / (Nfeeds3 - 1) # m
    Dl1 = D1 / wavelen
    Dl2 = D2 / wavelen
    Dl3 = D3 / wavelen
    for fd1 in range(Nfeeds1):
        cyl1_fp.append(fd1 * Dl1)
    for fd2 in range(Nfeeds2):
        cyl2_fp.append(fd2 * Dl2)
    for fd3 in range(Nfeeds3):
        cyl3_fp.append(fd3 * Dl3)
# Case 5: 3x32 D = 0.4 alternate 1/3
elif args.case == 5:
    Nfeeds = 32
    D = 0.4 # m
    Dl = D / wavelen
    for fd in range(Nfeeds):
        cyl1_fp.append((fd - 1.0/3) * Dl)
        cyl2_fp.append(fd * Dl)
        cyl3_fp.append((fd + 1.0/3) * Dl)
# Case 6: 3x32 10x0.8 + 1.0 + 20x1.1
elif args.case == 6:
    Nfeeds = 32
    D1 = 0.8 # m
    D2 = 1.0
    D3 = 1.1
    Dl1 = D1 / wavelen
    Dl2 = D2 / wavelen
    Dl3 = D3 / wavelen
    cyl1_sp = [Dl1] * 10 + [Dl3] * 5 + [Dl2] + [Dl3] * 15 # spacing between adjacent feeds
    cyl2_sp = [Dl3] * 10 + [Dl1] * 5 + [Dl2] + [Dl1] * 5 + [Dl3] * 10
    cyl3_sp = [Dl3] * 15 + [Dl2] + [Dl3] * 5 + [Dl1] * 10
    cyl1_fp = [fp for fp in np.cumsum(np.insert(cyl1_sp, 0, 0))]
    cyl2_fp = [fp for fp in np.cumsum(np.insert(cyl2_sp, 0, 0))]
    cyl3_fp = [fp for fp in np.cumsum(np.insert(cyl3_sp, 0, 0))]
else:
    raise Exception('Unsupported case: %d' % args.case)

cyls_fp = cyl1_fp + cyl2_fp + cyl3_fp # all feed v positions
Nf = len(cyls_fp) # total number of feeds
if args.auto_corr:
    v = [abs(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(i, Nf)]
else:
    v = [abs(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(i + 1, Nf)]
v = np.array(v)
v = np.around(v, 8)
v_uniq, v_cnt = unique(v, return_counts=True)

if args.outfile is not None:
    out_file = args.outfile
else:
    out_file = (conf[args.case -1] + '_%.1f_%s.hdf5') % (args.freq, args.auto_corr)

with h5py.File(out_file, 'w') as f:
    f.create_dataset('v', data=v_uniq)
    f.create_dataset('N', data=v_cnt)