#!/usr/bin/env python

"""Visualize uv-coverage.

:Authors: Shifan Zuo
:Date: 2014-08-01
:email: sfzuo@bao.ac.cn
:usage:
    python view_uv.py [-h] [-i IFREQ] [-a] [-c] [-o [OUTFILE]] [-f FIGFMT] [-l FIGLENGTH] [-w FIGWIDTH] [-g] [infile]
"""

import argparse

c = 3.0 * 10**8 # m/s, light speed

def visualize_uv(args):
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.infile, 'r') as f:
        baselines = f['baselines'][...]
        freqs = f['frequencies'][...]

    baselines = baselines.reshape(-1)
    if not args.all_freqs:
        freqs = np.array([freqs[args.ifreq]])

    bls, nus = np.broadcast_arrays(baselines[np.newaxis, :], freqs[:, np.newaxis])
    uv = bls * nus * 1.0 * 10**6 / c
    uv = uv.reshape(-1, 2)
    if args.conj:
        uv = np.concatenate((uv, -uv), axis=0)

    if not args.outfile:
        if args.all_freqs:
            outfile = 'uv_coverage_all_freqs.%s' % args.figfmt
        else:
            outfile = 'uv_coverage_%.1f.%s' % (freqs[args.ifreq], args.figfmt)

    # Plot uv-coverage
    plt.figure(figsize=(args.figlength, args.figwidth))
    plt.scatter(uv[:, 0], uv[:, 1])
    if args.grid:
        plt.grid()
    plt.xlabel('$u$')
    plt.ylabel('$v$')
    plt.savefig(outfile)



parser = argparse.ArgumentParser(description='Visualize uv-coverage.')
parser.add_argument('infile', type=str, nargs='?', help='Input hdf5 data file.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Single frequency channel index.')
parser.add_argument('-a', '--all_freqs', action='store_true', help='Plot all frequencies uv-coverage if true.')
parser.add_argument('-c', '--conj', action='store_false', help='Also include conjugate uv points.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into.')
parser.add_argument('-f', '--figfmt', default='pdf', help='Output image format.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add grids if true.')
parser.set_defaults(func=visualize_uv)

args = parser.parse_args()
args.func(args)

