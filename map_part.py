#!/usr/bin/env python

"""Get a square part of a full sky map.

:Authors: Shifan Zuo
:Date: 2014-08-26
:email: sfzuo@bao.ac.cn
:usage:
    python map_part.py [-h] [-o [OUTFILE]] [--clat [CLAT]] [--clon [CLON]] [--lat_ext LAT_EXT] [--lon_ext LON_EXT] [-i IFREQ] [-p {0,1,2,3}] inmap
"""

import argparse


def cart_proj(args):
    """Cartesian projection to get a square part of a full sky map.
    """
    import numpy as np
    import healpy
    import h5py
    import rotate as rot
    import matplotlib
    matplotlib.use('Agg')
    # from matplotlib import pyplot as plt

    # Read in maps data
    hpmap = None
    for mapname in args.mapfiles:
        with h5py.File(mapname, 'r') as f:
            if hpmap == None:
                hpmap = f['map'][:]
            else:
                hpmap += f['map'][:]

    lat_range = [args.clat - args.lat_ext, args.clat + args.lat_ext]
    lon_range = [args.clon - args.lon_ext, args.clon + args.lon_ext]
    # first rotate the map to let the point [clat, clon] be at the center
    # must rotate clon and clat separately, and NOTE the negative sign
    hpmap = rot.rotate_map(hpmap, rot=(-args.clon, 0.0, 0.0))
    hpmap = rot.rotate_map(hpmap, rot=(0.0, -args.clat, 0.0))

    # lat, lon value after rotation
    rot_lat_range = [lat_range[0] - args.clat, lat_range[1] - args.clat]
    rot_lon_range = [lon_range[0] - args.clon, lon_range[1] - args.clon]

    # cartesian projection
    cart_map = healpy.cartview(hpmap[args.ifreq, args.pol], latra=rot_lat_range, lonra=rot_lon_range, return_projected_map=True) # only T map

    out_file = args.outfile or 'cart_%.1f_%.1f_%.1f_%.1f.hdf5' % (lat_range[0], lat_range[1], lon_range[0], lon_range[1])
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=cart_map)
        f.attrs['lat_min'] = lat_range[0]
        f.attrs['lat_max'] = lat_range[1]
        f.attrs['lon_min'] = lon_range[0]
        f.attrs['lon_max'] = lon_range[1]


parser = argparse.ArgumentParser(description='Cartesian projection to get a square part of a full sky map.')
parser.add_argument('mapfiles', type=str, nargs='+', help='Input hdf5 sky map files to visualize. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the input args.')
parser.add_argument('--clat', type=float, nargs='?', default=-30, help='Central point latitude of the projected map.')
parser.add_argument('--clon', type=float, nargs='?', default=-90, help='Central point longitude of the projected map.')
parser.add_argument('--lat_ext', type=float, default=15.0, help='Latitude extension in either side of the central point.')
parser.add_argument('--lon_ext', type=float, default=15.0, help='Longitude extension in either side of the central point.')
parser.add_argument('-i', '--ifreq', type=int, default=0, help='Frequency channel index.')
parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.set_defaults(func=cart_proj)

args = parser.parse_args()
args.func(args)