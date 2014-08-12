#!/usr/bin/env python

"""Plot the :math:`a_{lm}` distribute of a healpix map.

:Authors: Shifan Zuo
:Date: 2014-06-16
:email: sfzuo@bao.ac.cn
:usage:
    python alm_dist.py [-h] [-o [OUT_FILE]] [-f FIGFMT] [-i [IFREQ]] [-p [POL]] [--maxl [MAXL]] [--lmin [LMIN]] [--lmax [LMAX]] [--mmin [MMIN]] [--mmax [MMAX]] [-l FIGLENGTH] [-w FIGWIDTH] [in_map]
"""

import argparse


def plot_alm(args):
    import os
    import numpy as np
    import h5py
    from cora.util import hputil
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.in_map, 'r') as f:
        in_map = f['map'][...]
    alm = hputil.sphtrans_sky(in_map, lmax=args.maxl)

    # plot alm
    plt.figure(figsize=(args.figlength, args.figwidth))
    plt.subplot(121)
    plt.pcolor(alm[args.ifreq, args.pol].T.real)
    lmin = args.lmin or 0
    lmax = args.lmax or alm.shape[-2]
    mmin = args.mmin or 0
    mmax = args.mmax  or alm.shape[-1]
    plt.xlim(lmin, lmax)
    plt.ylim(mmin, mmax)
    plt.xlabel(r'$l$')
    plt.ylabel(r'$m$')
    plt.title(r'$\Re\left(a_{lm}\right)$')
    plt.colorbar()

    plt.subplot(122)
    plt.pcolor(alm[args.ifreq, args.pol].T.imag)
    plt.xlim(lmin, lmax)
    plt.ylim(mmin, mmax)
    plt.xlabel(r'$l$')
    plt.ylabel(r'$m$')
    plt.title(r'$\Im\left(a_{lm}\right)$')
    plt.colorbar()
    out_file = args.out_file or ('alm_of_' + os.path.basename(args.in_map).replace('.hdf5', '_%d_%s.%s' % (args.ifreq, ('{%d}' % args.pol).format('T', 'Q', 'U', 'V'), args.figfmt)))
    plt.savefig(out_file)


parser = argparse.ArgumentParser(description='Plot the alm distribute of a healpix map.')
parser.add_argument('in_map', type=str, nargs='?', help='Input hdf5 sky map. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Name of the healpix map file to save.')
parser.add_argument('-f', '--figfmt', default='pdf', help='Output image format.')
parser.add_argument('-i', '--ifreq', type=int, nargs='?', default=0, help='Frequency channel index.')
parser.add_argument('-p', '--pol', type=int, nargs='?', default=0, help='Polarization index.')
parser.add_argument('--maxl', type=int, nargs='?', help='Max l in the spherical transform.')
parser.add_argument('--lmin', type=int, nargs='?', help='Min l to plot.')
parser.add_argument('--lmax', type=int, nargs='?', help='Max l to plot.')
parser.add_argument('--mmin', type=int, nargs='?', help='Min m to plot.')
parser.add_argument('--mmax', type=int, nargs='?', help='Max m to plot.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.set_defaults(func=plot_alm)

args = parser.parse_args()
args.func(args)
