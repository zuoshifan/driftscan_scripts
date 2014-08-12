#!/usr/bin/env python

"""Cutting a point source map at a given brightness temperature.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python pointsource_cut.py [-h] [-o [OUTFILE]] [-t THRESHOLD] [pointmap]
"""

import argparse


def cutting(args):
    """Cutting a point source map at a given brightness temperature.
    """
    import numpy as np
    import h5py

    # Read in map data
    with h5py.File(args.pointmap, 'r') as f:
        ptmap = f['map'][...]

    if args.threshold > 0:
        cut_map = np.where(ptmap<args.threshold, 0, ptmap)
    else:
        idx = np.unravel_index(np.argmax(ptmap), ptmap.shape) # the index of the max element
        cut_map = np.zeros_like(ptmap)
        cut_map[idx] = ptmap[idx]

    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    elif args.threshold > 0:
        out_file = ((args.pointmap.split('/')[-1]).split('.')[0]).replace('sim', 'cut') + '_' + str(int(args.threshold)) + '.hdf5'
    else:
        out_file = ((args.pointmap.split('/')[-1]).split('.')[0]).replace('sim', 'cut') + '_max.hdf5'

    # Save cut data
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=cut_map)

    print 'done!'



parser = argparse.ArgumentParser(description='Cut a point source map at a specified brightness temperature, zero is set below the cut threshold.')
parser.add_argument('pointmap', type=str, nargs='?', help='Input hdf5 point source sky map file.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the output hdf5 format healPix map.')
# parser.add_argument('-p', '--pol', type=int, default=0, choices=range(4), help='Polarization component to visualize, 0 for I/T, 1 for Q, 2 for U, 3 for V.')
parser.add_argument('-t', '--threshold', type=float, default=1000000.0, help='Cutting shreshold, below which will be set to zero. If less equal than 0, cutting threshold will be set to the max value of the input map.')
parser.set_defaults(func=cutting)

args = parser.parse_args()
args.func(args)