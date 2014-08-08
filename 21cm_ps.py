#!/usr/bin/env python

"""Compute neutral hydrogen 21cm temperature power spectrum.

:Authors: Shifan Zuo
:Date: 2014-04-08
:email: sfzuo@bao.ac.cn
:usage:
    python 
"""

import argparse
import numpy as np
import h5py

from cora.signal import corr21cm



def compute_ps(args):
    cr = corr21cm.Corr21cm()
    cr.ps_2d = False
    # 1d
    if not args.ps2d:
        if args.log:
            k_bands = np.logspace(np.log10(args.k_start), np.log10(args.k_end), args.k_num + 1, endpoint=True)
            k_type = 'log'
        else:
            k_bands = np.linspace(args.k_start, args.k_end, args.k_num +1, endpoint=True)
            k_type = 'lin'
        k_start = k_bands[:-1]
        k_end = k_bands[1:]
        k_center = 0.5 * (k_start + k_end)
        band_power = cr.ps_vv(k_center)
        delta2 = k_center**3 * band_power / (2 * np.pi**2)
        out_file = args.out_file or '21cm_ps1d.hdf5'
        with h5py.File(out_file, 'w') as f:
            f.create_dataset('k_bands', data=k_bands)
            f.create_dataset('k_start', data=k_start)
            f.create_dataset('k_end', data=k_end)
            f.create_dataset('k_center', data=k_center)
            f.create_dataset('powerspectrum', data=band_power)
            f.create_dataset('delta2', data=delta2)
            f.attrs['k_num'] = args.k_num
            f.attrs['k_type'] = k_type
    # 2d
    else:
        if args.log:
            kpar_bands = np.logspace(np.log10(args.kpar_start), np.log10(args.kpar_end), args.kpar_num + 1, endpoint=True)
            kperp_bands = np.logspace(np.log10(args.kperp_start), np.log10(args.kperp_end), args.kperp_num + 1, endpoint=True)
            k_type = 'log'
        else:
            kpar_bands = np.linspace(args.kpar_start, args.kpar_end, args.kpar_num + 1, endpoint=True)
            kperp_bands = np.linspace(args.kperp_start, args.kperp_end, args.kperp_num + 1, endpoint=True)
            k_type = 'lin'
        kpar_start = kpar_bands[:-1]
        kpar_end = kpar_bands[1:]
        kpar_center = 0.5 * (kpar_start + kpar_end)
        kperp_start = kperp_bands[:-1]
        kperp_end = kperp_bands[1:]
        kperp_center = 0.5 * (kperp_start + kperp_end)
        kpar_c, kperp_c = np.broadcast_arrays(kpar_center[np.newaxis, :], kperp_center[:, np.newaxis])
        kpar_c = kpar_c.flatten()
        kperp_c = kperp_c.flatten()
        k_center = (kpar_c**2 + kperp_c**2)**0.5
        band_power = cr.ps_vv(k_center)
        delta2 = kperp_c**2 * kpar_c * band_power / (2 * np.pi**2)
        out_file = args.out_file or '21cm_ps2d.hdf5'
        with h5py.File(out_file, 'w') as f:
            f.create_dataset('kpar_bands', data=kpar_bands)
            f.create_dataset('kpar_start', data=kpar_start)
            f.create_dataset('kpar_end', data=kpar_end)
            f.create_dataset('kperp_bands', data=kperp_bands)
            f.create_dataset('kperp_start', data=kperp_start)
            f.create_dataset('kperp_end', data=kperp_end)
            f.create_dataset('k_center', data=k_center)
            f.create_dataset('powerspectrum', data=band_power)
            f.create_dataset('delta2', data=delta2)
            f.attrs['kpar_num'] = args.kpar_num
            f.attrs['kperp_num'] = args.kperp_num
            f.attrs['k_type'] = k_type



parser = argparse.ArgumentParser(description='Compute neutral hydrogen 21cm temperature power spectrum.')
# 1d
parser.add_argument('--k_start', type=float, default=1.0e-3, help='Start k value for 1d power spectrum computation.')
parser.add_argument('--k_end', type=float, default=10.0, help='End k value for 1d power spectrum computation.')
parser.add_argument('--k_num', type=int, default=400, help='Number of k values for 1d power spectrum computation.')
# 2d
parser.add_argument('--kpar_start', type=float, default=0.01, help='Start k_parallel value for 2d power spectrum computation.')
parser.add_argument('--kpar_end', type=float, default=1.0, help='End k_parallel value for 2d power spectrum computation.')
parser.add_argument('--kpar_num', type=int, default=100, help='Number of k_parallel values for 2d power spectrum computation.')
parser.add_argument('--kperp_start', type=float, default=0.01, help='Start k_perpendicular value for 2d power spectrum computation.')
parser.add_argument('--kperp_end', type=float, default=1.0, help='End k_perpendicular value for 2d power spectrum computation.')
parser.add_argument('--kperp_num', type=int, default=100, help='Number of k_perpendicular values for 2d power spectrum computation.')
parser.add_argument('--ps2d', action='store_true', help='Compute 2d power spectrum if true.')
parser.add_argument('-o', '--out_file', help='Output power spectrum data file.')
parser.add_argument('--log', action="store_false", help='Use log scale if true.')
parser.set_defaults(func=compute_ps)

args = parser.parse_args()
args.func(args)
