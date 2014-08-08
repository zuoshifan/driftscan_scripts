#!/usr/bin/env python

"""Visualize a slice of beam transfer matrix :math:`(\mathbf{B}_{m}^{\nu})_{\alpha (Xl)}`.

:Authors: Shifan Zuo
:Date: 2014-05-26
:email: sfzuo@bao.ac.cn
:usage:
    python beam_matrix.py [-h] [-r [ROOT_DIR]] [-f [FILENAME]] [-m [MI]] [-o [OUT_DIR]] [-l FIGLENGTH] [-w FIGWIDTH]
"""

import argparse


def visualize_beam(args):
    import os
    from os.path import join
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    root_dir = args.root_dir if args.root_dir.endswith('/') else args.root_dir + '/'
    mi_len = len(str(args.mi))
    abs_filename = None
    for dirpath, subdirnames, filenames in os.walk(args.root_dir):
        for filename in filenames:
            if filename == args.filename and int(dirpath.split('/')[-1]) == args.mi:
                abs_filename = join(dirpath, filename)
                mi_subdir = dirpath.split('/')[-1]

    if abs_filename is None:
        raise Exception('No file %s for mi = %d' % (args.filename, args.mi))

    with h5py.File(abs_filename, 'r') as f:
        beam_m = f['beam_m'][...]
        freqs = f.attrs['frequencies']
    shape = beam_m.shape
    nfreq = shape[0]
    ntel = shape[1]*shape[2]
    nsky = shape[3]*shape[4]
    beam_m.shape = (nfreq, ntel, nsky) # reshape to :math:`(\mathbf{B}_{m}^{\nu})_{\alpha (Xl)}
    outdir = args.out_dir
    if not os.path.exists(outdir + mi_subdir):
        os.makedirs(outdir + mi_subdir)

    # plot
    for fi in range(nfreq):
        print '%d of %d...' % (fi, nfreq)
        plt.figure(figsize=(args.figlength, args.figwidth))
        plt.subplot(121)
        plt.pcolor(beam_m[fi].T.real)
        plt.xlim(0, ntel)
        plt.ylim(0, nsky)
        plt.xlabel(r"$\alpha$")
        plt.ylabel('$(Xl)$')
        plt.title(r'$\Re\left(\mathbf{B}_{\alpha (Xl)}\right)$ for $m = %d$ and $\nu = %.1f$ MHz' % (args.mi, freqs[fi]))
        plt.colorbar()
        
        plt.subplot(122)
        plt.pcolor(beam_m[fi].T.imag)
        plt.xlim(0, ntel)
        plt.ylim(0, nsky)
        plt.xlabel(r"$\alpha$")
        plt.ylabel('$(Xl)$')
        plt.title(r'$\Im\left(\mathbf{B}_{\alpha (Xl)}\right)$ for $m = %d$ and $\nu = %.1f$ MHz' % (args.mi, freqs[fi]))
        plt.colorbar()
        plt.savefig(outdir + mi_subdir + '/beam_%d_%d.png' % (args.mi, fi))


parser = argparse.ArgumentParser(description='Visualize a slice of beam transfer matrix.')
parser.add_argument('-r', '--root_dir', type=str, nargs='?', default='./', help='Input beam transfer matrix data files root directory.')
parser.add_argument('-f', '--filename', nargs='?', default='beam.hdf5', help='Input beam transfer matrix data file name.')
parser.add_argument('-m', '--mi', type=int, nargs='?', default=0, help='m index.')
parser.add_argument('-o', '--out_dir', nargs='?', default='beam/', help='Directory of output figures.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.set_defaults(func=visualize_beam)

args = parser.parse_args()
args.func(args)


