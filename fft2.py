#!/usr/bin/env python

"""2 dimensional FFT of a cartesian map.

:Authors: Shifan Zuo
:Date: 2014-08-26
:email: sfzuo@bao.ac.cn
:usage:
    python fft2.py [-h] [-o [OUTFILE]] inmap
"""

import argparse


def cart_fft2(args):
    """2 dimensional FFT of a cartesian map.
    """
    import numpy as np
    import h5py

    with h5py.File(args.inmap, 'r') as f:
        hpmap = f['map'][...]
    fft_map = np.fft.fft2(hpmap)
    out_file = args.outfile or 'fft2_' + args.inmap
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('fft', data=fft_map)


parser = argparse.ArgumentParser(description='2 dimensional FFT of a cartesian map.')
parser.add_argument('inmap', type=str, help='Input healpix map.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the input args.')
parser.set_defaults(func=cart_fft2)

args = parser.parse_args()
args.func(args)