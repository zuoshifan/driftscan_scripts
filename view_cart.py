#!/usr/bin/env python

"""Visualize a cartesian sky map.

:Authors: Shifan Zuo
:Date: 2014-08-26
:email: sfzuo@bao.ac.cn
:usage:
    python view_cart.py [-h] [-o [OUTFILE]] [--figfmt FIGFMT] [--vmin VMIN] [--vmax VMAX] [-l FIGLENGTH] [-w FIGWIDTH] [-g] inmap
"""

import argparse


def visualize_cart(args):
    """Visualize a cartesian sky map.
    """
    import numpy as np
    import healpy
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator

    with h5py.File(args.inmap, 'r') as f:
        cart_map = f['map'][...]
        lat_min = f.attrs['lat_min']
        lat_max = f.attrs['lat_max']
        lon_min = f.attrs['lon_min']
        lon_max = f.attrs['lon_max']

    if cart_map.ndim > 2:
        cart_map = cart_map[args.ifreq, args.pol]
        out_file = args.outfile or ((args.inmap.split('/')[-1]).split('.')[0] + '_' + str(args.ifreq) + '_{' + str(args.pol) + '}' +  '.' + args.figfmt).format('T', 'Q', 'U', 'V')
    else:
        out_file = args.outfile or args.inmap.split('/')[-1].replace('.hdf5', '.' + args.figfmt)

    plt.figure(figsize=(args.figlength, args.figwidth))
    # ext = (lon_max + 360, lon_min + 360, lat_min, lat_max)
    ext = (lon_max, lon_min, lat_min, lat_max)
    plt.imshow(cart_map, extent=ext, origin='lower', vmin=args.vmin, vmax=args.vmax)
    # cart_map = cart_map[:, 315-15:315+15]
    # plt.imshow(cart_map, extent=ext, aspect='auto', origin='lower', vmin=args.vmin, vmax=args.vmax)
    # plt.imshow(cart_map, origin='lower', vmin=args.vmin, vmax=args.vmax)
    plt.xlabel('RA / deg')
    plt.ylabel('DEC / deg')
    cbar = plt.colorbar()
    if args.figfmt == 'pdf':
        cbar.solids.set_rasterized(True) # eliminate lines in colorbar for PDF output, see http://matplotlib.1069221.n5.nabble.com/rasterized-colorbar-td39582.html
    ## can not completely eliminate lines in colorbar, see http://matplotlib.org/api/pyplot_api.html?highlight=colorbar#matplotlib.pyplot.colorbar
    # cbar.solids.set_edgecolor("face")
    # plt.draw()
    ## -----------------------------------------------------
    plt.savefig(out_file, bbox_inches='tight')


parser = argparse.ArgumentParser(description='Visualize a cartesian sky map.')
parser.add_argument('inmap', type=str, help='Input cartesian map.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the input args.')
parser.add_argument('--figfmt', default='png', help='Output image format.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('--vmin', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--vmax', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=visualize_cart)

args = parser.parse_args()
args.func(args)
