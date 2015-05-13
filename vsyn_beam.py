#!/usr/bin/env python

"""Synthesized beam in `v` direction computation.

:Authors: Shifan Zuo
:Date: 2015-05-12
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
parser.add_argument('-c', '--case', type=int, choices=range(1, 7), default=1, help='Which array configuration.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-o', '--outfile', help='Output name file name.')
args = parser.parse_args()




conf = ['3x32_0.5', '3x32_1.0', '3x32_2.5', '31_32_33_1.0', '31_32_33_2.5', 'random_1.0']
Nconf = len(conf)
cyl1_fp = [] # feed v position for feed 1, unit: lambda
cyl2_fp = []
cyl3_fp = []

# Case 1: 3x32 v = 0.5
if args.case == 1:
    Nfeeds = 32
    v = 0.5 # in unit of one wavelength
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * v)
        cyl2_fp.append(fd * v)
        cyl3_fp.append(fd * v)
# Case 2: 3x32 v = 1.0
elif args.case == 2:
    Nfeeds = 32
    v = 1.0 # in unit of one wavelength
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * v)
        cyl2_fp.append(fd * v)
        cyl3_fp.append(fd * v)
# Case 3: 3x32 v = 2.5
elif args.case == 3:
    Nfeeds = 32
    v = 2.5 # in unit of one wavelength
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * v)
        cyl2_fp.append(fd * v)
        cyl3_fp.append(fd * v)
# Case 4: 31+32+33 v = 1.0
elif args.case == 4:
    Nfeeds1 = 31
    Nfeeds2 = 32
    Nfeeds3 = 33
    v2 = 1.0 # in unit of one wavelength
    v1 = (Nfeeds2 - 1) * v2 / (Nfeeds1 - 1)
    v3 = (Nfeeds2 - 1) * v2 / (Nfeeds3 - 1)
    for fd1 in range(Nfeeds1):
        cyl1_fp.append(fd1 * v1)
    for fd2 in range(Nfeeds2):
        cyl2_fp.append(fd2 * v2)
    for fd3 in range(Nfeeds3):
        cyl3_fp.append(fd3 * v3)
# Case 5: 31+32+33 v = 2.5
elif args.case == 5:
    Nfeeds1 = 31
    Nfeeds2 = 32
    Nfeeds3 = 33
    v2 = 2.5 # in unit of one wavelength
    v1 = (Nfeeds2 - 1) * v2 / (Nfeeds1 - 1)
    v3 = (Nfeeds2 - 1) * v2 / (Nfeeds3 - 1)
    for fd1 in range(Nfeeds1):
        cyl1_fp.append(fd1 * v1)
    for fd2 in range(Nfeeds2):
        cyl2_fp.append(fd2 * v2)
    for fd3 in range(Nfeeds3):
        cyl3_fp.append(fd3 * v3)
# Case 6: random v = 1.0
elif args.case == 6:
    Nfeeds = 32
    v = 1.0 # in unit of one wavelength
    for fd in range(Nfeeds):
        cyl1_fp.append(fd * v + np.random.uniform(-0.5, 0.5))
        cyl2_fp.append(fd * v + np.random.uniform(-0.5, 0.5))
        cyl3_fp.append(fd * v + np.random.uniform(-0.5, 0.5))



cyls_fp = cyl1_fp + cyl2_fp + cyl3_fp # all feed v positions
Nf = len(cyls_fp) # total number of feeds
if args.auto_corr:
    v_list = [(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(Nf)]
else:
    v_list = [(cyls_fp[i] - cyls_fp[j]) for i in range(Nf) for j in range(Nf) if i != j]
# v = np.array(v)

def vsynbeam(theta):
    syn = 0.0
    for v in v_list:
        syn += np.exp(2 * np.pi * 1.0J * v * np.sin(np.radians(theta)))

    return syn.real / len(v_list)

if args.outfile is not None:
    out_file = args.outfile
else:
    out_file = (conf[args.case -1] + '_vsyn_beam_%s.%s') % (args.auto_corr, args.figfmt)

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
