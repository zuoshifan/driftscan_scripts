#!/usr/bin/env python

"""Synthesized beam in `v` direction computation only for regular array.

:Authors: Shifan Zuo
:Date: 2015-05-13
:email: sfzuo@bao.ac.cn
"""


import argparse
import numpy as np
# import healpy
# import h5py
# from mpi4py import MPI
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator


# Read in arguments
parser = argparse.ArgumentParser(description="Calculate v baselines distribution.")
parser.add_argument('-c', '--case', type=int, choices=range(1, 4), default=1, help='Which array configuration.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-o', '--outfile', help='Output name file name.')
args = parser.parse_args()




conf = ['3x32_0.5', '3x32_1.0', '3x32_2.5']
Nconf = len(conf)
cyl1_fp = [] # feed v position for feed 1, unit: lambda
cyl2_fp = []
cyl3_fp = []

# Case 1: 3x32 v = 0.5
if args.case == 1:
    Nfeeds = 32
    v = 0.5 # in unit of one wavelength
# Case 2: 3x32 v = 1.0
if args.case == 2:
    Nfeeds = 32
    v = 1.0 # in unit of one wavelength
# Case 3: 3x32 v = 2.5
if args.case == 3:
    Nfeeds = 32
    v = 2.5 # in unit of one wavelength


# cyls_fp = cyl1_fp + cyl2_fp + cyl3_fp # all feed v positions
# Nf = len(cyls_fp) # total number of feeds
# if args.auto_corr:
#     v_list = [(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(Nf)]
# else:
#     v_list = [(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(Nf) if i != j]
# # v = np.array(v)

def vsynbeam(theta):
    syn = 0.0
    for n in range(Nfeeds):
        syn += 2 * (Nfeeds - n) * np.cos(2 * np.pi * n * v * np.sin(np.radians(theta)))

    return syn / Nfeeds**2

if args.outfile is not None:
    out_file = args.outfile
else:
    out_file = (conf[args.case -1] + '_vsyn_beam_Alt.%s') % (args.figfmt)

theta = np.linspace(-90.0, 90.0, 10001)
plt.figure(figsize=(args.figlength,args.figwidth))
plt.plot(theta, vsynbeam(theta))
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
plt.xlabel(r'$\theta$ / degree')
plt.xlim(-90.0, 90.0)
plt.ylim(0.0, 1.0)
plt.grid()
plt.savefig(out_file)
