#!/usr/bin/env python

"""Visualize `lm` slice of beam transfer matrix.

:Authors: Shifan Zuo
:Date: 2014-04-13
:email: sfzuo@bao.ac.cn
:usage:
    python beam_lm.py [-h] [-f [FILENAME]] [-i IFREQ] [-p {0,1,2,3}] [-b BL BL] [--mmax MMAX] [--lmax LMAX] [-a] [-l FIGLENGTH] [-w FIGWIDTH] [indir]
"""

import argparse


def visualize_beam_lm(args):
    import os
    from os.path import join
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    beam_lms = []
    mi_idx = []
    for dirpath, subdirnames, filenames in os.walk(args.indir):
        for filename in filenames:
            if filename == args.filename:
                mi = int(dirpath.split('/')[-1])
                mi_idx.append(mi)
                abs_filename = join(dirpath, filename)
                with h5py.File(abs_filename, 'r') as f:
                    beam_m = f['beam_m'][...]
                    # print 'original shape: ', beam_m.shape
                    freq = f.attrs['frequencies'][args.ifreq]
                    if args.single_bl:
                        beam_m = beam_m[args.ifreq, args.bl[0], args.bl[1], args.pol, :]
                        # print 'slice shape: ', beam_m.shape
                    else:
                        beam_m = np.sum(np.sum(beam_m, axis=1), axis=1)[args.ifreq, args.pol, :]
                        # print 'summed shape: ', beam_m.shape
                    beam_lms.append(beam_m)
    beam_lms = np.array(beam_lms)
    beam_lms.take(mi_idx, axis=0)  # sort according to the `m` index
    beam_lms = beam_lms[0:args.mmax, 0:args.lmax]

    # Plot the image
    plt.figure(figsize=(args.figlength, args.figwidth))
    plt.subplot(121)
    # plt.imshow(beam_lms.T.real, origin='lower')  # abscissa m, ordinate l
    plt.imshow(beam_lms.real, origin='lower')  # abscissa l, ordinate m
    plt.xlabel('$l$')
    plt.ylabel('$m$')
    plt.title('$B_{l,m}$ real')
    plt.colorbar()
    plt.subplot(122)
    # plt.imshow(beam_lms.T.imag, origin='lower')  # abscissa m, ordinate l
    plt.imshow(beam_lms.imag, origin='lower')  # abscissa l, ordinate m
    plt.xlabel('$l$')
    plt.ylabel('$m$')
    plt.title('$B_{l,m}$ image')
    plt.colorbar()
    if args.single_bl:
        plt.savefig('beam_lm_%d_%d_%d_%d_%d_%d.png' % (freq, args.pol, args.bl[0], args.bl[1], beam_lms.shape[1], beam_lms.shape[0]))  # lmax, mmax
    else:
        plt.savefig('beam_lm_%d_%d_all_bl_%d_%d.png' % (freq, args.pol, beam_lms.shape[1], beam_lms.shape[0]))  # lmax, mmax


parser = argparse.ArgumentParser(description='Visualize `lm` slice of beam transfer matrix.')
parser.add_argument('indir', type=str, nargs='?', help='Input beam transfer matrix data files directory.')
parser.add_argument('-f', '--filename', nargs='?', default='beam.hdf5', help='Input beam transfer matrix data file name.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('-b', '--bl', type=int, nargs=2, default=[0, 0], help='Baseline index.')
parser.add_argument('--mmax', type=int, help='Max m to plot.')
parser.add_argument('--lmax', type=int, help='Max l to plot.')
parser.add_argument('-a', '--single_bl', action='store_true', help='Plot single baseline `lm` slice of beam transfer matrix if present, else the sum of all baselines.')
# parser.add_argument('--log', type=bool, action='store_true', help='Plot log10(B_lm) when present.')
# parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
# parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.set_defaults(func=visualize_beam_lm)

args = parser.parse_args()
args.func(args)

