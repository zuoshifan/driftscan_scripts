#!/usr/bin/env python

"""Scatter plot a sky map in hdf5 files.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python scatter_view.py [-h] [-o [OUTFILE]] [-f FIGFMT] [-i IFREQ] [-p {0,1,2,3}] [--min MIN] [--max MAX] [-w FIGWIDTH] [-t FIGHEIGHT] [-g] mapfiles [mapfiles ...]
"""

import argparse


def scatter_plt(args):
    """Scatter plot a set of data.
    """
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    # Read in maps data
    hpmap = None
    for mapname in args.mapfiles:
        with h5py.File(mapname, 'r') as f:
            if hpmap == None:
                hpmap = f['map'][:]
            else:
                hpmap += f['map'][:]

    # Check args validity
    if args.ifreq < 0 and args.ifreq >= (hpmap.shape)[0]:
        raise Exception('Invalid frequency channel %d, should be in range(0, %d).'%(args.ifreq, (hpmap.shape)[0]))
    if args.pol >= (hpmap.shape)[1]:
        raise Exception('Invalid frequency channel %d, should be in range(0, %d).'%(args.pol, (hpmap.shape)[1]))
    if args.figwidth <= 0:
        raise Exception('Image width figwidth (= %f) must greater than 0'%args.figwidth)
    if args.figheight <= 0:
        raise Exception('Image height figheight (= %f) must greater than 0'%args.figheight)

    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    else:
        out_file = ((args.mapfiles[0].split('/')[-1]).split('.')[0] + '_' + str(args.ifreq) + '_{' + str(args.pol) + '}' +  '.' + args.figfrm).format('T', 'Q', 'U', 'V')

    # Create a scatter plot of the data
    mean = np.mean(hpmap[args.ifreq][args.pol])
    std = np.std(hpmap[args.ifreq][args.pol])
    plt.figure(figsize=(args.figwidth, args.figheight))
    plt.plot(np.arange((hpmap.shape)[-1]), hpmap[args.ifreq][args.pol], 'o')
    plt.title('mean = %f, std = %f'%(mean, std))
    plt.ylim(args.min, args.max)
    plt.ylabel('K')
    if args.grid:
        plt.grid()
    plt.savefig('sct_' + out_file)


parser = argparse.ArgumentParser(description='Scatter plot a sky map in hdf5 files.')
parser.add_argument('mapfiles', type=str, nargs='+', help='Input hdf5 sky map files to visualize. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto create from the first input map file name and input args.')
parser.add_argument('-f', '--figfmt', default='pdf', help='Output image format.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0).')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-w', '--figwidth', type=float, default=8, help='Output image width.')
parser.add_argument('-t', '--figheight', type=float, default=6, help='Output image height.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=scatter_plt)

args = parser.parse_args()
args.func(args)
