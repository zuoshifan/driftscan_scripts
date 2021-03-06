#!/usr/bin/env python

"""Visualize the cylinder array beam.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python cylinder_beam.py [-h] [-o [OUTFILE]] [--lat [LAT]] [--lon [LON]] [-f [FREQ]] [-d [CYL_WIDTH]] [-n [NSIDE]] [--feed [FEED]] [--min MIN] [--max MAX] [-l FIGLENGTH] [-w FIGWIDTH] [-g]
"""

import argparse


def visualize_cyl_beam(args):
    """Visualize cylinder beam.
    """
    import numpy as np
    import healpy
    import hpvisual
    from drift.telescope import cylinder as cyl
    # from drift.core import telescope as tel
    from cora.util import hputil
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator

    if args.pol:
        pol_cyl = cyl.PolarisedCylinderTelescope(args.lat, args.lon)
        # pol_cyl.zenith = tel.latlon_to_sphpol([args.lat, args.lon]) # NOTE: this is wrong, use the constructor above or the following assignment, because the assignment will automatically trigger drift/util/config.py __set__() method, which will call tel.latlon_to_sphpol.
        # pol_cyl.zenith = [args.lat, args.lon]
        pol_cyl.cylinder_width = args.cyl_width
        pol_cyl._frequencies = np.array([args.freq])
        pol_cyl._angpos = hputil.ang_positions(args.nside)
        beamx_resp = pol_cyl.beamx(args.feed, 0)
        beamy_resp = pol_cyl.beamy(args.feed, 0)
        # Create output image file name
        if args.outfile:
            out_file = args.outfile
        else:
            freq = ('%s'%args.freq).replace('.', '_')
            cyl_width = ('%s'%args.cyl_width).replace('.', '_')
            out_file = 'pol_cyl_beam_resp_%s__%s.%s'%(freq, cyl_width, args.figfmt)
        # Plot and save image
        fig1 = plt.figure(1, figsize=(args.figlength, 2*args.figwidth))
        title11 = 'X feed theta direction, f = %s MHz, w = %s m'%(args.freq, args.cyl_width)
        title12 = 'X feed phi direction, f = %s MHz, w = %s m'%(args.freq, args.cyl_width)
        hpvisual.mollview(beamx_resp[:, 0], fig=1, sub=221, title=title11, min=args.min, max=args.max)
        hpvisual.mollview(beamx_resp[:, 1], fig=1, sub=222, title=title12, min=args.min, max=args.max)
        title21 = 'Y feed theta direction, f = %s MHz, w = %s m'%(args.freq, args.cyl_width)
        title22 = 'Y feed phi direction, f = %s MHz, w = %s m'%(args.freq, args.cyl_width)
        hpvisual.mollview(beamy_resp[:, 0], fig=1, sub=223, title=title21, min=args.min, max=args.max)
        hpvisual.mollview(beamy_resp[:, 1], fig=1, sub=224, title=title22, min=args.min, max=args.max)
        if args.grid:
            healpy.graticule()
        fig1.savefig(out_file)
        fig1.clf()

        # R_{I to I}
        II_resp = (beamx_resp[:, 0]**2 + beamx_resp[:, 1]**2 + beamy_resp[:, 0]**2 + beamy_resp[:, 1]**2) / 2.0
        # R_{P to I}
        PI_resp = (beamx_resp[:, 0]**2 - beamx_resp[:, 1]**2 + beamx_resp[:, 0]*beamx_resp[:, 1] + beamy_resp[:, 0]**2 - beamy_resp[:, 1]**2 + beamy_resp[:, 0]*beamy_resp[:, 1]) / 2.0

        # Plot the R_{I to I}response
        fig2 = plt.figure(2, figsize=(args.figlength, args.figwidth))
        hpvisual.mollview(II_resp, fig=2, min=args.min, max=args.max)
        if args.grid:
            healpy.graticule()
        fig2.savefig('II_mollview' + out_file)
        fig2.clf()

        # Plot the R_{P to I}response
        fig3 = plt.figure(3, figsize=(args.figlength, args.figwidth))
        hpvisual.mollview(PI_resp, fig=3, min=args.min, max=args.max)
        if args.grid:
            healpy.graticule()
        fig3.savefig('PI_mollview' + out_file)
        fig3.clf()


        # cart_array = healpy.cartview(II_resp, fig=3, lonra=[-10, 10], min=args.min, max=args.max, aspect=0.25, notext=True, hold=True, title='', return_projected_map=True)
        lon_range = [-10, 10]
        ext = (-10, 10, -90, 90)
        # R_{I to I } cartesian
        II_cart = hpvisual.cartview(II_resp, lonra=lon_range, return_projected_map=True)
        ticks = np.linspace(np.min(II_cart), np.max(II_cart), 6)
        ticks = np.around(ticks, decimals=1)
        plt.figure(figsize=(4, 5))
        plt.imshow(II_cart, extent=ext, aspect=0.2, origin='lower', vmin=ticks[0], vmax=ticks[-1])
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(30))
        plt.xlabel('EW / deg')
        plt.ylabel('NS / deg')
        cbar = plt.colorbar(shrink=0.6 ,orientation='horizontal', ticks=ticks)
        cbar.solids.set_rasterized(True)
        plt.savefig('II_cartview_' + out_file)
        # R_{P to I } cartesian
        PI_cart = hpvisual.cartview(PI_resp, lonra=lon_range, return_projected_map=True)
        ticks = np.linspace(np.min(PI_cart), np.max(PI_cart), 4)
        ticks = np.around(ticks, decimals=2)
        plt.figure(figsize=(4, 5))
        plt.imshow(PI_cart, extent=ext, aspect=0.2, origin='lower', vmin=ticks[0], vmax=ticks[-1])
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(30))
        plt.xlabel('EW / deg')
        plt.ylabel('NS / deg')
        cbar = plt.colorbar(shrink=0.6 ,orientation='horizontal', ticks=ticks)
        cbar.solids.set_rasterized(True)
        plt.savefig('PI_cartview_' + out_file)

    else:
        unpol_cyl = cyl.UnpolarisedCylinderTelescope(args.lat, args.lon)
        unpol_cyl.cylinder_width = args.cyl_width
        unpol_cyl._frequencies = np.array([args.freq])
        unpol_cyl._angpos = hputil.ang_positions(args.nside)
        beam_amp = unpol_cyl.beam(args.feed, 0)
        # Create output image file name
        if args.outfile:
            out_file = args.outfile
        else:
            freq = ('%s'%args.freq).replace('.', '_')
            cyl_width = ('%s'%args.cyl_width).replace('.', '_')
            out_file = 'unpol_cyl_beam_resp_%s__%s.%s'%(freq, cyl_width, args.figfmt)
        # Plot and save image
        fig = plt.figure(1, figsize=(args.figlength, args.figwidth))
        title = 'Beam pattern, f = %s MHz, w = %s m'%(args.freq, args.cyl_width)
        hpvisual.mollview(beam_amp, fig=1, sub=111, coord='G', title=title, min=args.min, max=args.max)
        if args.grid:
            healpy.graticule()
        fig.savefig(out_file)
        fig.clf()



parser = argparse.ArgumentParser(description='Visualize cylinder telescope beam response.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the input args.')
parser.add_argument('--figfmt', default='png', help='Output image format.')
parser.add_argument('--lat', type=float, nargs='?', default=45, help='Telescope latitude.')
parser.add_argument('--lon', type=float, nargs='?', default=0, help='Telescope longitude.')
parser.add_argument('-f', '--freq', type=float, nargs='?', default=700.0, help='Telescope observing frequency.')
parser.add_argument('-d', '--cyl_width', type=float, nargs='?', default=15.0, help='Cylinder width.')
parser.add_argument('-n', '--nside', type=int, nargs='?', default=256, help='Healpix nside.')
parser.add_argument('-p', '--pol', action='store_false', help='Polarized beam response if present, else unpolarized.')
parser.add_argument('--feed', type=int, nargs='?', default=0, help='Which feed.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=visualize_cyl_beam)

args = parser.parse_args()
args.func(args)