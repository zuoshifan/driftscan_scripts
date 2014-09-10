#!/usr/bin/env python

"""Plot the :math:`a_{lm}` distribute of a healpix map for plane wave expansion.

:Authors: Shifan Zuo
:Date: 2014-09-09
:email: sfzuo@bao.ac.cn
:usage:
    python palm_dist.py [-h] [-o [OUT_FILE]] [-f FIGFMT] [--lmin [LMIN]] [--lmax [LMAX]] [--mmin [MMIN]] [--mmax [MMAX]] [--vmin VMIN] [--vmax VMAX] [--log] [-l FIGLENGTH] [-w FIGWIDTH] [in_file]
"""

import argparse


def plot_alm(args):
    import os
    import numpy as np
    import numpy.ma as ma
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.in_file, 'r') as f:
        alm = f['alm'][...]

    # mask |m| > l values
    mask = np.zeros_like(alm, dtype=np.bool)
    for row in range(mask.shape[0]):
        for col in range(mask.shape[1]):
            if abs(col - (alm.shape[1] - 1)/2) > row:
                mask[row, col] = 1
    mask.dtype = bool
    alm = ma.array(alm, mask=mask) # mask array

    # plot alm
    lmin = args.lmin or 0
    lmax = args.lmax or alm.shape[-2]
    mmin = args.mmin or - (alm.shape[-1] - 1)/2
    mmax = args.mmax  or (alm.shape[-1] - 1)/2
    ext = (lmin, lmax, mmin, mmax)
    plt.figure(figsize=(args.figlength, args.figwidth))
    if args.log:
        plt.imshow(np.log(np.abs(alm.T)), extent=ext, aspect='auto', origin='lower', vmin=args.vmin, vmax=args.vmax)
    else:
        plt.imshow(np.abs(alm.T), extent=ext, aspect='auto', origin='lower', vmin=args.vmin, vmax=args.vmax)
    plt.xlabel(r'$l$')
    plt.ylabel(r'$m$')
    # plt.title(r'$\Re\left(a_{lm}\right)$')
    cbar = plt.colorbar()
    if args.figfmt == 'pdf':
        cbar.solids.set_rasterized(True)

    out_file = args.out_file or args.in_file.replace('hdf5', args.figfmt)
    plt.savefig(out_file)


parser = argparse.ArgumentParser(description='Plot the alm distribute of a healpix map.')
parser.add_argument('in_file', type=str, nargs='?', help='Input hdf5 spherical coefficients file.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Name of the healpix map file to save.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
parser.add_argument('--lmin', type=int, nargs='?', help='Min l to plot.')
parser.add_argument('--lmax', type=int, nargs='?', help='Max l to plot.')
parser.add_argument('--mmin', type=int, nargs='?', help='Min m to plot.')
parser.add_argument('--mmax', type=int, nargs='?', help='Max m to plot.')
parser.add_argument('--vmin', type=float, help='Min value to plot.')
parser.add_argument('--vmax', type=float, help='Max value to plot.')
parser.add_argument('--log', action='store_true', help='Whether plot log scale.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.set_defaults(func=plot_alm)

args = parser.parse_args()
args.func(args)








