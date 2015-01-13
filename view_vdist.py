#!/usr/bin/env python

"""Visualize v baselines distribution.

:Authors: Shifan Zuo
:Date: 2015-01-08
:email: sfzuo@bao.ac.cn
"""

import argparse
import numpy as np


# def natpattern(n):
#     """Pattern that prints out a number upto `n` (natural number - no sign)."""
#     return ("%0" + repr(int(np.ceil(np.log10(n + 1)))) + "d")

def visualize_vdist(args):
    # import os
    # from os.path import join
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt


    with h5py.File(args.data_file, 'r') as f:
        v = f['v'][...]
        N = f['N'][...]

    # plot
    plt.figure(figsize=(args.figlength, args.figwidth))
    plt.plot(v, N, 'o')
    # plt.xlim(0.0, 12.5)
    if args.min and args.max:
        plt.ylim(args.min, args.max)
    plt.xlabel(r'$v\ /\ \lambda$')
    plt.ylabel(r'$N$')
    plt.grid()
    if args.outfile is not None:
        outfile = args.outfile
    else:
        outfile = args.data_file.replace('.hdf5', '.' + args.figfmt)
    plt.savefig(outfile)


parser = argparse.ArgumentParser(description='Visualize visibilities.')
parser.add_argument('data_file', type=str, nargs='?', help='Input data file.')
parser.add_argument('--figfmt', default='png', help='Output image format.')
# parser.add_argument('--log', type=bool, action='store_true', help='Plot log10(B_lm) when present.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-o', '--outfile', help='Output name file name.')
parser.set_defaults(func=visualize_vdist)

args = parser.parse_args()
args.func(args)
