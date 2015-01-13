#!/usr/bin/env python

"""Visualize a single frequency sky temperature (no polarization) map in hdf5 files.

:Authors: Shifan Zuo
:Date: 2014-12-06
:email: sfzuo@bao.ac.cn
"""

import argparse


def visualize_map(args):
    """Visualize a single frequency sky temperature (no polarization) map in hdf5 files.

    Arguments
    ---------
    args : argparse namespace.
    """
    import numpy as np
    import h5py
    import healpy
    import hpvisual
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    # Read in map data
    hpmap = None
    for mapname in args.mapfiles:
        with h5py.File(mapname, 'r') as f:
            if hpmap == None:
                hpmap = f['map'][:]
            else:
                hpmap += f['map'][:]

    # Check args validity
    if args.figlength <= 0:
        raise Exception('Figure length figlength (= %f) must greater than 0'%args.figlength)
    if args.figwidth <= 0:
        raise Exception('Figure width figwidth (= %f) must greater than 0'%args.figwidth)

    # chose appropriate data
    if args.data == 'a':
        hpmap = np.abs(hpmap)
        data = 'abs'
    elif args.data == 'r':
        hpmap = hpmap.real
        data = 'real'
    else:
        hpmap = hpmap.imag
        data = 'imag'
    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    else:
        out_file = mapname.replace('.hdf5', '_'+data+'.'+args.figfmt)

    # Plot and save image
    fig = plt.figure(1, figsize=(args.figlength,args.figwidth))
    hpvisual.mollview(hpmap, fig=1, title='', min=args.min, max=args.max)
    # healpy.cartview(hpmap[ifreq][args.pol], fig=1, title='', min=args.min, max=args.max)
    # fig = plt.figure()
    # ax = fig.add_axes()
    # cbar.solids.set_rasterized(True)
    if args.grid:
        healpy.graticule()
    if args.tight:
        fig.savefig(out_file, bbox_inches='tight')
    else:
        fig.savefig(out_file)
    fig.clf()


parser = argparse.ArgumentParser(description='Visualize a sky map in hdf5 files.')
parser.add_argument('mapfiles', type=str, nargs='+', help='Input hdf5 sky map files to visualize. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the first input map file name and input args.')
parser.add_argument('-d', '--data', type=str, choices=['a', 'r', 'i'], default='a', help='Which data to visualize, `a` for abs, `r` for real and `i` for image.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
# parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
# parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.add_argument('-t', '--tight', action='store_true', help='Tight the figure marin space.')
parser.set_defaults(func=visualize_map)

args = parser.parse_args()
args.func(args)