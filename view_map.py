#!/usr/bin/env python

"""Visualize sky maps in hdf5 files.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python view_map.py [-h] [-o [OUTFILE]] [-f FIGFMT] [-i IFREQ] [-p {0,1,2,3}] [--min MIN] [--max MAX] [-l FIGLENGTH] [-w FIGWIDTH] [-g] mapfiles [mapfiles ...]
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
    import hpvisual
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
        out_file = ((args.mapfiles[0].split('/')[-1]).split('.')[0] + '_' + str(ifreq) + '_{' + str(args.pol) + '}' +  '.' + args.figfmt).format('T', 'Q', 'U', 'V')

    # Plot and save image
    if args.view == 'o':
        fig = plt.figure(1, figsize=(8, 6))
    else:
        fig = plt.figure(1, figsize=(args.figlength,args.figwidth))

    map_data = hpmap[ifreq][args.pol]
    if args.sqrt:
        map_data = map_data / np.sqrt(np.abs(map_data))
        map_data = map_data / np.sqrt(np.abs(map_data))
        # map_data = map_data / np.sqrt(np.abs(map_data))

    # smoothing the map with a Gaussian symmetric beam
    if args.fwhm is not None:
        fwhm = np.radians(args.fwhm)
        map_data = healpy.smoothing(map_data, fwhm=fwhm)

    if args.view == 'm':
        # set color map
        if args.cmap is None:
            cmap = None
        else:
            if args.cmap == 'jet09':
                import colormap
                cmap = colormap.jet09
            else:
                from pylab import cm
                # cmap = cm.hot
                cmap = getattr(cm, args.cmap)
            cmap.set_under('w')

        if args.abs:
            healpy.mollview(np.abs(map_data), fig=1, title='', cmap=cmap, min=args.min, max=args.max)
        else:
            healpy.mollview(map_data, fig=1, title='', cmap=cmap, min=args.min, max=args.max)
    elif args.view == 'c':
        healpy.cartview(map_data, fig=1, title='', min=args.min, max=args.max)
    elif args.view == 'o':
        hpvisual.orthview(map_data, rot=(0, 90, 0), fig=1, title='', min=args.min, max=args.max, half_sky=True) # rot to make NCP at the center
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
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
parser.add_argument('-v', '--view', type=str, choices=['m', 'c', 'o'], default='m', help='Which view, `m` for mollview, `c` for cartview, `o` for orthview.')
parser.add_argument('-s', '--sqrt', action='store_true', help='Plot by sqrt of the map.')
parser.add_argument('--fwhm', type=float, default=None, help='Smoothing the map with a Gaussian symmetric beam with FWHM this value, in degree.')
parser.add_argument('-a', '--abs', action='store_true', help='Plot the abs value of the map.')
parser.add_argument('-c', '--cmap', type=str, default=None, help='The cmap to use.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel to visualize (start from 0). Negative integer N means the last Nth channel.')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.add_argument('-t', '--tight', action='store_true', help='Tight the figure marin space.')
parser.set_defaults(func=visualize_map)

args = parser.parse_args()
args.func(args)
