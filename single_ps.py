#!/usr/bin/env python

"""Make a single pointsource healpy map at a specified position (theta, phi).

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python single_ps.py [-h] [-o [OUTFILE]] [--nside NSIDE] [--nfreq NFREQ] [--freq_lower FREQ_LOWER] [--freq_upper FREQ_UPPER] [--T T] [--npol {1,4}] [theta] [phi]
"""

import argparse


def mk_single_ps(args):
    """Make a single pointsource healpy map.
    
    Arguments:
    - `args`:
    """
    import numpy as np
    import h5py
    import healpy
    
    ps_map = np.zeros((args.nfreq, args.npol, 12*args.nside**2))
    pix = healpy.ang2pix(args.nside, np.radians(args.theta), np.radians(args.phi))
    for ifreq in range(args.nfreq):
        ps_map[ifreq, 0, pix] = args.T

    out_file = args.outfile or 'single_ps_%d_%d_%d_%d_%d_%d.hdf5'%(args.nside, args.freq_lower, args.freq_upper, args.nfreq, args.theta, args.phi)
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=ps_map)


parser = argparse.ArgumentParser(description='Make a single pointsource healpy map at a specified position (theta, phi).')
parser.add_argument('theta', type=float, nargs='?', help='Theta angular coordinates (in degree) of a point on the sphere.')
parser.add_argument('phi', type=float, nargs='?', help='Phi angular coordinates (in degree) of a point on the sphere.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the output image file (png/eps) to save into. If not present, the output image file name will be auto created from the input args (png).')
parser.add_argument('--nside', type=int, default=256, help='Healpix map nside.')
parser.add_argument('--nfreq', type=int, default=5, help='Number of frequencies channel.')
parser.add_argument('--freq_lower', type=float, default=700.0, help='Lower freqency.')
parser.add_argument('--freq_upper', type=float, default=800.0, help='Upper freqency.')
parser.add_argument('--T', type=float, default=10000.0, help='Brightness temperature of the point source.')
parser.add_argument('--npol', type=int, default=4, choices=[1,4], help='Polarization 4 or unpolarization 1.')
parser.set_defaults(func=mk_single_ps)

args = parser.parse_args()
args.func(args)