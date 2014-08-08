#!/usr/bin/env python

"""Plot the SVD spectrum of the :math:`\bar{B}` matrix.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python svd_spectrum.py [-h] [-o OUTFILE] [-i IFREQ] [-l FIGLENGTH] [-w FIGWIDTH] [-n AUTOLEVELS] [--lmin LMIN] [--lmax LMAX] [--lint LINT] [filename]
"""

import argparse


def plt_svd(args):
    import h5py
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.filename, 'r') as f:
        svs = f['singularvalues'][...]
        svs_norm = svs / np.max(svs)
        
    # Check args validity
    if args.ifreq < -(svs.shape)[1] or args.ifreq >= (svs.shape)[1]:
        raise Exception('Invalid frequency channel %d, should be in range(-%d, %d).'%(args.ifreq, (svs.shape)[1], (svs.shape)[1]))
    else:
        ifreq = args.ifreq if args.ifreq >= 0 else args.ifreq + (svs.shape)[1]
    if args.figlength <= 0:
        raise Exception('Figure length figlength (= %f) must greater than 0'%args.figlength)
    if args.figwidth <= 0:
        raise Exception('Figure width figwidth (= %f) must greater than 0'%args.figwidth)

    # plot and save figure
    plt.figure(figsize=(args.figlength, args.figwidth))
    cs1 = plt.contourf(np.log10(svs_norm[:, ifreq, :]).T, args.autolevels)
    lmin = args.lmin or np.ceil(np.min(cs1.levels))
    lmax = args.lmax or np.floor(np.max(cs1.levels))
    levels = [lmin] if lmin == lmax else np.arange(lmin, lmax, args.lint)
    cs2 = plt.contour(cs1, levels=levels, colors='k', linestyles='solid', hold='on')
    plt.clabel(cs2, inline=1)
    cbar = plt.colorbar(cs1)
    cbar.add_lines(cs2)
    plt.title('SVD spectrum - $\log_{10}(\Sigma / \Sigma_{max})$')
    plt.xlabel('m')
    plt.ylabel('SVD mode')
    outfile = args.outfile or 'svd_spectrum_%d.png'%ifreq
    plt.savefig(outfile)


parser = argparse.ArgumentParser(description='Plot the SVD spectrum of the B bar matrix.')
parser.add_argument('filename', type=str, nargs='?', default='svdspectrum.hdf5', help='Input hdf5 svdspectrum file.')
parser.add_argument('-o', '--outfile', help='Output image file name.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Out figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-n', '--autolevels', type=int, default=200, help='contour N automatically-chosen levels.')
parser.add_argument('--lmin', type=float, help='Min level curves to draw.')
parser.add_argument('--lmax', type=float, help='Max level curves to draw.')
parser.add_argument('--lint', type=float, default=2, help='Interval of level curves.')
parser.set_defaults(func=plt_svd)

args = parser.parse_args()
args.func(args)