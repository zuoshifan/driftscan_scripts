#!/usr/bin/env python

"""Visualize a dish telescope beam response.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python view_beam.py [-h] [-o [OUTFILE]] [--figfmt FIGFMT] [--lat [LAT]] [--lon [LON]] [-f [FREQ]] [-d [DISH_WIDTH]] [-n [NSIDE]] [--feed [FEED]] [--min MIN] [--max MAX] [-l FIGLENGTH] [-w FIGWIDTH] [-g]
"""

import argparse
import numpy as np
import healpy
from scipy.special import jn
from cora.util import coord
from cora.util import hputil, units
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt



def wavelengths(frequencies):
    """The central wavelength of each frequency band (in metres)."""
    return units.c / (1e6 * frequencies)

def latlon_to_sphpol(latlon):
    zenith = np.array([np.pi / 2.0 - np.radians(latlon[0]),
                       np.remainder(np.radians(latlon[1]), 2*np.pi)])
    return zenith

def beam_circular(angpos, zenith, diameter):
    """Beam pattern for a uniformly illuminated circular dish.

    Parameters
    ----------
    angpos : np.ndarray
        Array of angular positions
    zenith : np.ndarray
        Co-ordinates of the zenith.
    diameter : scalar
        Diameter of the dish (in units of wavelength).

    Returns
    -------
    beam : np.ndarray
        Beam pattern at each position in angpos.
    """

    def jinc(x):
        return 0.5 * (jn(0, x) + jn(2, x))

    x = (1.0 - coord.sph_dot(angpos, zenith)**2)**0.5 * np.pi * diameter

    return 2*jinc(x)


# Implement the X and Y beam patterns (assuming all feeds are identical).
# These need to return a vector for each position on the sky
# (self._angpos) in thetahat, phihat coordinates.
def beamx(feed, freq, angpos, zenith, dish_width):
    # Calculate beam amplitude
    beam = beam_circular(angpos, zenith, dish_width / wavelengths(freq))

    # Add a vector direction to beam - X beam is EW (phihat)
    beam = beam[:, np.newaxis] * np.array([0.0, 1.0])

    return beam

def beamy(feed, freq, angpos, zenith, dish_width):
    # Calculate beam amplitude
    beam = beam_circular(angpos, zenith, dish_width / wavelengths(freq))

    # Add a vector direction to beam - Y beam is NS (thetahat)
    # Fine provided beam does not cross a pole.
    beam = beam[:, np.newaxis] * np.array([1.0, 0.0])

    return beam



def visualize_beam(args):
    angpos = hputil.ang_positions(args.nside)
    zenith = latlon_to_sphpol([args.lat, args.lon])
    beamx_resp = beamx(args.feed, args.freq, angpos, zenith, args.dish_width)
    beamy_resp = beamy(args.feed, args.freq, angpos, zenith, args.dish_width)
    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    else:
        freq = ('%s'%args.freq).replace('.', '_')
        diameter = ('%s'%args.dish_width).replace('.', '_')
        out_file = 'beam_resp_%s__%s.%s'%(freq, diameter, args.figfmt)
    # Plot and save image
    fig = plt.figure(1, figsize=(args.figlength, args.figwidth))
    title1 = 'X beam pattern, f= %sMHz, d = %sm'%(args.freq, args.dish_width)
    healpy.mollview(beamx_resp[:, 1], fig=1, sub=121, title=title1, min=args.min, max=args.max)
    title2 = 'Y beam pattern, f= %sMHz, d = %sM'%(args.freq, args.dish_width)
    healpy.mollview(beamy_resp[:, 0], fig=1, sub=122, title=title2, min=args.min, max=args.max)
    if args.grid:
        healpy.graticule()
    fig.savefig(out_file)
    fig.clf()



parser = argparse.ArgumentParser(description='Visualize a dish telescope beam response.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the input args.')
parser.add_argument('--figfmt', default='pdf', help='Output image format.')
parser.add_argument('--lat', type=float, nargs='?', default=30, help='Telescope latitude.')
parser.add_argument('--lon', type=float, nargs='?', default=0, help='Telescope longitude.')
parser.add_argument('-f', '--freq', type=float, nargs='?', default=100.0, help='Telescope observing frequency.')
parser.add_argument('-d', '--dish_width', type=float, nargs='?', default=3.5, help='Dish diameter.')
parser.add_argument('-n', '--nside', type=int, nargs='?', default=256, help='Healpix nside.')
parser.add_argument('--feed', type=int, nargs='?', default=0, help='Which feed.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=visualize_beam)

args = parser.parse_args()
args.func(args)