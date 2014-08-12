#!/usr/bin/env python

"""Plot the S/N or S/F spectrum of the KL transform.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python snr_spectrum.py [-h] [-o OUTFILE] [-l FIGLENGTH] [-w FIGWIDTH] [-n AUTOLEVELS] [--lminsn LMINSN] [--lmaxsn LMAXSN] [--lintsn LINTSN] [--lminsf LMINSF] [--lmaxsf LMAXSF] [--lintsf LINTSF] [-k KLMODES] [filename]
"""

import argparse


def plt_sn(args):
    import h5py
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.filename, 'r') as f:
        evals = f['evals'][...]
        f_evals = None # only double KL has f_evals
        try:
            f_evals = f['f_evals'][...]
        except:
            pass
    if f_evals == None:
        plt.figure(figsize=(args.figlength, args.figwidth))
        cs1 = plt.contourf(np.log10(evals).T[::-1][:], args.autolevels)
        lminsn = args.lminsn or np.ceil(np.min(cs1.levels))
        lmaxsn = args.lmaxsn or np.floor(np.max(cs1.levels))
        levels = [lminsn] if lminsn == lmaxsn else np.arange(lminsn, lmaxsn, args.lintsn)
        cs2 = plt.contour(cs1, levels=levels, colors='k', linestyles='solid', hold='on')
        plt.clabel(cs2, inline=1)
        cbar = plt.colorbar(cs1)
        cbar.add_lines(cs2)
        plt.title('Signal/Noise Ratio - log10(S/N)')
        plt.xlabel('m')
        plt.ylabel('KL mode number (sorted)')
        if args.klmodes:
            plt.ylim(0, args.klmodes)
        outfile = args.outfile or 'kl_evals.%s' % args.figfmt
        plt.savefig(outfile)
    else:
        plt.figure(figsize=(2 * args.figlength, args.figwidth))
        plt.subplot(121)
        cs1 = plt.contourf(np.log10(f_evals).T[::-1][:], args.autolevels)
        lminsf = args.lminsf or np.ceil(np.min(cs1.levels))
        lmaxsf = args.lmaxsf or np.floor(np.max(cs1.levels))
        levels = [lminsf] if lminsf == lmaxsf else np.arange(lminsf, lmaxsf, args.lintsf)
        cs2 = plt.contour(cs1, levels=levels, colors='k', linestyles='solid', hold='on')
        plt.clabel(cs2, inline=1)
        cbar = plt.colorbar(cs1)
        cbar.add_lines(cs2)
        plt.title('Signal/Foreground Ratio - log10(S/F)')
        plt.xlabel('m')
        plt.ylabel('KL mode number (sorted)')
        if args.klmodes:
            plt.ylim(0, args.klmodes)
        plt.subplot(122)
        cs3 = plt.contourf(np.log10(evals).T[::-1][:], args.autolevels)
        lminsn = args.lminsn or np.ceil(np.min(cs3.levels))
        lmaxsn = args.lmaxsn or np.floor(np.max(cs3.levels))
        levels = [lminsn] if lminsn == lmaxsn else np.arange(lminsn, lmaxsn, args.lintsn)
        cs4 = plt.contour(cs3, levels=levels, colors='k', linestyles='solid', hold='on')
        plt.clabel(cs4, inline=1)
        cbar = plt.colorbar(cs3)
        cbar.add_lines(cs4)
        plt.title('Signal/Noise Ratio - log10(S/N)')
        plt.xlabel('m')
        plt.ylabel('KL mode number (sorted)')
        if args.klmodes:
            plt.ylim(0, args.klmodes)
        outfile = args.outfile or 'dk_evals.%s' % args.figfmt
        plt.savefig(outfile)



parser = argparse.ArgumentParser(description='Plot the S/N or S/F spectrum of the KL transform.')
parser.add_argument('filename', type=str, nargs='?', default='evals.hdf5', help='Input hdf5 evals file.')
parser.add_argument('-o', '--outfile', help='Output image file name.')
parser.add_argument('-f', '--figfmt', default='pdf', help='Output image format.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Out figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-n', '--autolevels', type=int, default=200, help='contour N automatically-chosen levels.')
parser.add_argument('--lminsn', type=float, help='Min level curves to draw in S/N evals.')
parser.add_argument('--lmaxsn', type=float, help='Max level curves to draw in S/N evals.')
parser.add_argument('--lintsn', type=float, default=2, help='Interval of level curves in S/N evals.')
parser.add_argument('--lminsf', type=float, help='Min level curves to draw in S/F evals.')
parser.add_argument('--lmaxsf', type=float, help='Max level curves to draw in S/F evals.')
parser.add_argument('--lintsf', type=float, default=2, help='Interval of level curves in S/F evals.')
parser.add_argument('-k', '--klmodes', type=int, help='Max kl modes to plot.')
parser.set_defaults(func=plt_sn)

args = parser.parse_args()
args.func(args)