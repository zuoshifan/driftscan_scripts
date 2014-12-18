#!/usr/bin/env python

"""Visualize visibilities waterfall plot.

:Authors: Shifan Zuo
:Date: 2014-12-08
:email: sfzuo@bao.ac.cn
"""

import argparse
import warnings
import numpy as np


def natpattern(n):
    """Pattern that prints out a number upto `n` (natural number - no sign)."""
    return ("%0" + repr(int(np.ceil(np.log10(n + 1)))) + "d")

def visualize_vis(args):
    import os
    from os.path import join
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt


    data_file = args.indir + '/ts_commondata.hdf5'
    with h5py.File(data_file, 'r') as f:
        baselines = f['baselines'][...]
        tphi = f['phi'][...]
        freqs = f['frequencies'][...]
    print 'Number of baselines: ', baselines.shape[0]
    u = baselines[args.bl, 0] # m, EW direction
    v = baselines[args.bl, 1] # m, NS direction
    phi = np.linspace(0, 2*np.pi, tphi.shape[0]+1, endpoint=True)

    nfreq = freqs.shape[0]
    delta_nu = freqs[1] - freqs[0]
    freq_lower = freqs[0] - 0.5 * delta_nu
    freq_upper = freqs[-1] + 0.5 * delta_nu
    # freqs = freq_lower + (np.arange(nfreq) + 0.5) * ((freq_upper - freq_lower) / nfreq)
    freqs_bound = np.linspace(freq_lower, freq_upper, nfreq + 1, endpoint=True)

    args.ifreq = args.ifreq if args.ifreq >= 0 else args.ifreq + nfreq
    ifreq_str = natpattern(nfreq) % args.ifreq

    ts = []
    ifreqs = []
    for dirpath, subdirnames, filenames in os.walk(args.indir):
        for filename in filenames:
            if filename.startswith('timestream_'):
                ifreq = int((filename.split('.')[0].split('_')[-1]))
                ifreqs.append(ifreq)
                abs_filename = join(dirpath, filename)
                with h5py.File(abs_filename, 'r') as f:
                    ts.append(f['timestream'][args.bl])
    isort = np.argsort(ifreqs)
    # ifreqs = np.take(ifreqs, isort)
    # sort ts for ascending frequencies
    ts = np.array(ts)
    ts1 = np.zeros_like(ts)
    cnt = 0
    for idx in isort:
        ts1[cnt] = ts[idx]
        cnt += 1

    # plot
    plt.figure(figsize=(args.figlength, args.figwidth))
    plt.subplot(211)
    plt.pcolor(freqs_bound, phi, ts.T[::-1, :].real)
    if args.pmin is None or args.pmin < phi[0]:
        pmin = phi[0]
    else:
        pmin = args.pmin
    if args.pmax is None or args.pmax > phi[-1]:
        pmax = phi[-1]
    else:
        pmax = args.pmax
    if pmin > pmax:
        warnings.warn('pmin > pmax, exchange the two value')
    if pmin == pmax:
        raise Exception('pmin == pmax')
    plt.ylim(pmax, pmin)
    plt.xlabel(r'$\nu$ / MHz')
    plt.ylabel(r'$\phi$ / rad')
    cbar = plt.colorbar()
    if args.figfmt == 'pdf':
        cbar.solids.set_rasterized(True)

    plt.subplot(212)
    plt.pcolor(freqs_bound, phi, ts.T[::-1, :].imag)
    plt.ylim(pmax, pmin)
    plt.xlabel(r'$\nu$ / MHz')
    plt.ylabel(r'$\phi$ / rad')
    cbar = plt.colorbar()
    if args.figfmt == 'pdf':
        cbar.solids.set_rasterized(True)

    plt.savefig('visibility_2D_%.1f_%.1f_%d_%.1f_%.1f.%s' % (u, v, args.bl, pmin, pmax, args.figfmt))



parser = argparse.ArgumentParser(description='Visualize visibilities.')
parser.add_argument('indir', type=str, nargs='?', help='Input visibility time stream data files directory.')
# parser.add_argument('-f', '--filename', nargs='?', default='beam.hdf5', help='Input beam transfer matrix data file name.')
parser.add_argument('--figfmt', default='png', help='Output image format.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-b', '--bl', type=int, nargs='?', default=0, help='Baseline index.')
parser.add_argument('--pmin', type=float, help='Lowest phi range.')
parser.add_argument('--pmax', type=float, help='Highest phi range.')
# parser.add_argument('--log', type=bool, action='store_true', help='Plot log10(B_lm) when present.')
# parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
# parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.set_defaults(func=visualize_vis)

args = parser.parse_args()
args.func(args)
