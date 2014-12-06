#!/usr/bin/env python

"""compute the spherical harmonic coefficients of plane wave in equatorial coordinate.

:Authors: Shifan Zuo
:Date: 2014-09-09
:email: sfzuo@bao.ac.cn
:usage:
    python planewave.py [-h] [-c {1,2}] [-a] [-f FREQ] [-n NSIDE] [-o OUT_FILE]
"""

import argparse
import numpy as np
import scipy.special as special
import healpy
import h5py

from cora.util import hputil


# Read in arguments
parser = argparse.ArgumentParser(description="compute the spherical harmonic coefficients of plane wave in equatorial coordinate.")
parser.add_argument('-c', '--case', type=int, choices=[1, 2], default=1, help='Which array configuration, 1 for 32+32+32 case, 2 for 31+32+33 case.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-f', '--freq', type=float, default=750.0, help='Observing frequency.')
# parser.add_argument('-n', '--nside', type=int, default=256, help='Healpix nside.')
parser.add_argument('-b', '--baseline', type=float, nargs=2, help='A particular baseline in unit m in format [u, v] in local ground coordinate. Use all baselines of an array if not given.')
parser.add_argument('-o', '--out_file', help='Output name file name.')
args = parser.parse_args()


def vec_len(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def vec_ang(vec):
    leng = vec_len(vec)
    if leng < 1.0e-8:
        return 0.0, 0.0

    vec /= leng # normalize
    x, y, z = vec
    theta = np.arccos(z)
    if np.abs(x) < 1.0e-8:
        phi = np.pi /2.0 if y >= 0.0 else -np.pi/2.0
    else:
        if x >= 0.0 and y >= 0.0:
            N = 0
        elif x < 0.0:
            N = 1
        else:
            N = 2
        phi = np.arctan(vec[1] / vec[0]) + N * np.pi
    return theta, phi

def sph_jn(n, z):
    """Spherical Bessel function of real int order n.

    See: http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html"""
    if z < 1.0e-12:
        return 1.0 if n == 0 else 0.0

    return np.sqrt(0.5 * np.pi / z) * special.jn(n+0.5, z)


case = 0

ncyl = 3
cyl_w = 15.0 # m

# observing frequency
freq = args.freq # MHz
c = 3.0e8 # m/s
wavelen = c / (freq * 1.0e6)

# transform matrix from local coordinate (`x` points toward East, `y` points toward North, `z` points toward the Zenith) to equatorial coordinate (`X` points toward the local meridian in the `yz` plane, `Z` points toward the North pole, `Y` is perpendicular to `X` and `Z` and they together construct a right hand coordinate)
R = np.array([[0.0, -np.sqrt(2)/2.0, np.sqrt(2)/2.0],
              [1.0, 0.0, 0.0],
              [np.sqrt(2)/2.0, 0.0, np.sqrt(2)/2.0]])


if args.baseline is None:
    # feeds positions in local coordinate
    local_pos = []
    # Case 1: 32 + 32 + 32
    if args.case == 1:
        case = args.case
        print 'Case%d' % case
        nfeeds = 32
        feeds_spacing = 0.4 # m
        for cyl in range(ncyl):
            for feed in range(nfeeds):
                local_pos.append(np.array([cyl * cyl_w, feed * feeds_spacing, 0.0]))
    # Case 2: 31 + 32 + 33
    elif args.case == 2:
        case = args.case
        print 'Case%d' % case
        nfeeds = [31, 32, 33]
        feeds_spacing = [12.4/30, 0.4, 12.4/32]
        for cyl in range(ncyl):
            for feed in range(nfeeds[cyl]):
                local_pos.append(np.array([cyl * cyl_w, feed * feeds_spacing[cyl], 0.0]))

    # feeds positions in equatorial coordinate
    eq_pos = [np.dot(R, pos) for pos in local_pos]
    pos_len = len(eq_pos)
    # baseline vectors in equatorial coordinate
    auto_corr = args.auto_corr
    if auto_corr:
        bls = [(eq_pos[pos1] - eq_pos[pos2]) for pos1 in range(pos_len) for pos2 in range(pos1, pos_len)]
    else:
        bls = [(eq_pos[pos1] - eq_pos[pos2]) for pos1 in range(pos_len) for pos2 in range(pos1+1, pos_len)]

        uvecs = [(bl / wavelen) for bl in bls] # lists of np.array([u, v, w])
    # ulen = [vec_len(u) for u in uvecs] # lists of |u|
    # uang = [vec_ang(u) for u in uvecs] # list of (theta, phi)
else:
    bls = [np.dot(R, np.array([args.baseline[0], args.baseline[1], 0.0]))] # transform to equatorial coordinate
    uvecs = [(bl / wavelen) for bl in bls] # lists of np.array([u, v, w])

N = 32 # nfeeds
D = 0.4 # m
lmax = 2 * np.pi * np.sqrt((ncyl * cyl_w)**2 + (N * D)**2) / wavelen
mmax = 2 * np.pi * ncyl * cyl_w / wavelen
lmax = np.int(np.ceil(lmax))
mmax = np.int(np.ceil(mmax))

alm = np.zeros((lmax+1, 2*mmax+1), dtype=np.complex128)
# sum over all baselines
for u in uvecs:
    ulen = vec_len(u)
    theta, phi = vec_ang(u)
    # print 'ulen = %f, theta = %f, phi = %f' % (ulen, theta, phi)
    for l in range(lmax+1):
        print '%d of %d...' % (l, lmax)
        jl = sph_jn(l, 2 * np.pi * ulen)
        # print 'jl = %f' % jl
        for m in range(-l, l+1):
            alm[l, m] += 4 * np.pi * jl * special.sph_harm(m, l, phi, theta)

# save data
outfile = args.out_file or 'palm_%.1f_%s_case%d.hdf5' % (args.freq, args.auto_corr, case)
with h5py.File(outfile, 'w') as f:
    f.create_dataset('alm', data=alm)



# # healpix
# nside = args.nside
# nx, ny, nz = healpy.pix2vec(nside, range(12*nside**2))
# enus = np.zeros_like(nx, dtype=np.complex128)
# for uvec in uvecs:
#     # enus += np.exp(2*np.pi*1.0J*np.abs(uvec[0] * nx + uvec[1] * ny + uvec[2] * nz))
#     enus += np.exp(2*np.pi*1.0J*(uvec[0] * nx + uvec[1] * ny + uvec[2] * nz))

# # spherical harmonic transform
# alm = hputil.sphtrans_complex(enus, centered=True)

# # save data
# outfile = args.out_file or 'palm_%d_%.1f_%s_case%d.hdf5' % (args.nside, args.freq, args.auto_corr, args.case)
# with h5py.File(outfile, 'w') as f:
#     f.create_dataset('alm', data=alm)
