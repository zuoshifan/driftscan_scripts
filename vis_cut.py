#!/usr/bin/env python

"""Cut a time slice of the visibilities.

:Authors: Shifan Zuo
:Date: 2015-04-05
:email: sfzuo@bao.ac.cn
"""

import argparse


def cut(args):
    """Cut a time slice of the visibilities.

    Arguments
    ---------
    args : argparse namespace.
    """
    import numpy as np
    import h5py

    # format the slice arg
    if len(args.slice) == 1:
        slc = [args.slice[0], args.slice[0]+1]
    elif len(args.slice) == 2:
        slc = args.slice
        print slc
    else:
        raise ValueError('Invalid slice value')

    # Read in maps data
    for vis_file in args.visfiles:
        with h5py.File(vis_file, 'r+') as f:
            timestream = f['timestream'][...]
            vis_slice = f['timestream'][:, slc[0]:slc[1]]
            f['timestream'][...] = np.zeros(timestream.shape, dtype=timestream.dtype)
            f['timestream'][:, slc[0]:slc[1]] = vis_slice


parser = argparse.ArgumentParser(description='Cut a time slice of the visibilities.')
parser.add_argument('visfiles', type=str, nargs='+', help='Input hdf5 visibilities files. If more than one, they will be processed one by one.')
parser.add_argument('-s', '--slice', nargs='*', type=int, default=[0], help='The time slice of the visibilities we will get (zeros will be set for others). If no args, use the default value; if one args given, get the give value slice; If two given, get the range between the two; error otherwise.')

parser.set_defaults(func=cut)

args = parser.parse_args()
args.func(args)
