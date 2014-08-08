#!/usr/bin/env python

"""Visualize GMRT feeds position configuration.

:Authors: Shifan Zuo
:Date: 2014-05-12
:email: sfzuo@bao.ac.cn
:usage:
    python gmrt_cong.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Read in data
feeds_pos_file = '/home/zuoshifan/programming/python/21cmcosmology/driftscan/drift/telescope/gmrtpositions.dat'
feeds_pos = np.loadtxt(feeds_pos_file)

# plot
plt.figure(figsize=(8, 6))
plt.scatter(feeds_pos[:, 0], feeds_pos[:, 1])
plt.xlabel('E')
plt.ylabel('N')
plt.savefig('gmrt_conf.png')

