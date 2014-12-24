#!/usr/bin/env python

"""Visualize matrix B^{dagger}B, :math:`\left(\mathbf{B}^{\dagger}\mathbf{B}\right)_{(Xl)(X'l')}`.

:Authors: Shifan Zuo
:Date: 2014-06-21
:email: sfzuo@bao.ac.cn
"""

import argparse


def natpattern(n):
    """Pattern that prints out a number upto `n` (natural number - no sign)."""
    return ("%0" + repr(int(np.ceil(np.log10(n + 1)))) + "d")


def visualize_bb(args):
    import os
    # from os.path import join
    import re
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    from drift.util import mpiutil

    # root_dir = args.root_dir if args.root_dir.endswith('/') else args.root_dir + '/'
    # mi_len = len(str(args.mi))
    # abs_filename = None
    # for dirpath, subdirnames, filenames in os.walk(args.root_dir):
    #     for filename in filenames:
    #         if filename == args.filename and int(dirpath.split('/')[-1]) == args.mi:
    #             abs_filename = join(dirpath, filename)
    #             mi_subdir = dirpath.split('/')[-1]

    # if abs_filename is None:
    #     raise Exception('No file %s for mi = %d' % (args.filename, args.mi))

    mi = re.search('[0-9]+', args.filename).group(0)
    mi_subdir = mi

    with h5py.File(args.filename, 'r') as f:
        beam_m = f['beam_m'][...]

    pardir = os.path.abspath(os.path.join(os.path.dirname(args.filename), os.pardir)) # get parent dir
    with h5py.File(pardir + '/telescope_data.hdf5', 'r') as f:
        freqs = f['frequencies'][...]

    shape = beam_m.shape
    nfreq = shape[0]
    ntel = shape[1]*shape[2]
    nsky = shape[3]*shape[4]
    beam_m.shape = (nfreq, ntel, nsky) # reshape to :math:`(\mathbf{B}_{m}^{\nu})_{\alpha (Xl)}
    bb_m = np.zeros((nfreq, nsky, nsky), dtype=beam_m.dtype)
    for fi in range(nfreq):
        bb_m[fi] = np.dot(beam_m[fi].T.conj(), beam_m[fi])
    outdir = args.out_dir
    if not os.path.exists(outdir + mi_subdir):
        os.makedirs(outdir + mi_subdir)

    # plot
    for fi in mpiutil.mpirange(nfreq):
        print '%d of %d...' % (fi, nfreq)
        plt.figure(figsize=(args.figlength, args.figwidth))
        plt.subplot(121)
        plt.pcolor(bb_m[fi].real, vmin=args.min, vmax=args.max)
        plt.xlim(0, nsky)
        plt.ylim(0, nsky)
        plt.xlabel(r"$(X'l')$")
        plt.ylabel('$(Xl)$')
        plt.title(r"$\Re\left[\left(\mathbf{B}^{\dagger}\mathbf{B}\right)_{(Xl) (X'l')}\right]$ for $m = %d$ and $\nu = %.1f$ MHz" % (int(mi), freqs[fi]))
        plt.colorbar()

        plt.subplot(122)
        plt.pcolor(bb_m[fi].imag, vmin=args.min, vmax=args.max)
        plt.xlim(0, nsky)
        plt.ylim(0, nsky)
        plt.xlabel(r"$(X'l')$")
        plt.ylabel('$(Xl)$')
        plt.title(r"$\Im\left[\left(\mathbf{B}^{\dagger}\mathbf{B}\right)_{(Xl) (X'l')}\right]$ for $m = %s$ and $\nu = %.1f$ MHz" % (mi, freqs[fi]))
        plt.colorbar()
        plt.savefig(outdir + mi_subdir + '/bb_%s_%d.%s' % (mi, fi, args.figfmt))
        plt.close()


parser = argparse.ArgumentParser(description='Visualize matrix B^{dagger}B.')
# parser.add_argument('-r', '--root_dir', type=str, nargs='?', default='./', help='Input beam transfer matrix data files root directory.')
parser.add_argument('filename', nargs='?', default='beam.hdf5', help='Input beam transfer matrix data file name.')
# parser.add_argument('-m', '--mi', type=int, nargs='?', default=0, help='m index.')
parser.add_argument('-o', '--out_dir', nargs='?', default='bb/', help='Directory of output figures.')
parser.add_argument('--figfmt', default='png', help='Output image format.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.set_defaults(func=visualize_bb)

args = parser.parse_args()
args.func(args)
