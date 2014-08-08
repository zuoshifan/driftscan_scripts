#!/usr/bin/env python

"""Visualize sky maps in hdf5 files.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python view_map.py [-h] [-o [OUTFILE]] [-i IFREQ] [-p {0,1,2,3}] [--min MIN] [--max MAX] [-l FIGLENGTH] [-w FIGWIDTH] [-g] mapfiles [mapfiles ...]
"""

import argparse


def visualize_map(args):
    """Visualize sky maps in hdf5 files.
    
    Arguments
    ---------
    args : argparse namespace.
    """
    import numpy as np
    import h5py
    import healpy
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
    if args.ifreq < -(hpmap.shape)[0] or args.ifreq >= (hpmap.shape)[0]:
        raise Exception('Invalid frequency channel %d, should be in range(-%d, %d).'%(args.ifreq, (hpmap.shape)[0], (hpmap.shape)[0]))
    else:
        ifreq = args.ifreq if args.ifreq >= 0 else args.ifreq + (hpmap.shape)[0]
    if args.pol >= (hpmap.shape)[1]:
        raise Exception('Invalid frequency channel %d, should be in range(0, %d).'%(args.pol, (hpmap.shape)[1]))
    if args.figlength <= 0:
        raise Exception('Figure length figlength (= %f) must greater than 0'%args.figlength)
    if args.figwidth <= 0:
        raise Exception('Figure width figwidth (= %f) must greater than 0'%args.figwidth)

    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    else:
        out_file = ((args.mapfiles[0].split('/')[-1]).split('.')[0] + '_' + str(ifreq) + '_{' + str(args.pol) + '}' +  '.png').format('T', 'Q', 'U', 'V')

    # Plot and save image
    fig = plt.figure(1, figsize=(args.figlength,args.figwidth))
    healpy.mollview(hpmap[ifreq][args.pol], fig=1, title='', min=args.min, max=args.max)
    # healpy.cartview(hpmap[ifreq][args.pol], fig=1, title='', min=args.min, max=args.max)
    if args.grid:
        healpy.graticule()
    fig.savefig(out_file)
    fig.clf()


parser = argparse.ArgumentParser(description='Visualize a sky map in hdf5 files.')
parser.add_argument('mapfiles', type=str, nargs='+', help='Input hdf5 sky map files to visualize. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file (png/eps) to save into. If not present, the output image file name will be auto created from the first input map file name and input args (png).')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=visualize_map)

args = parser.parse_args()
args.func(args)