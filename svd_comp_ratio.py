#!/usr/bin/env python

"""SVD compression ratio.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage:
    python svd_comp_ratio.py
"""

import os
from os.path import join, getsize

# root_dir = 'disharray2/'
root_dir = ''
mode_size = 0
svd_size = 0
for root, dirs, files in os.walk(root_dir + 'yamldriver/timestream/mmodes/'):
    mode_size += sum(getsize(join(root, name)) for name in files if name == 'mode.hdf5')
    svd_size += sum(getsize(join(root, name)) for name in files if name == 'svd.hdf5')
print 'Total size of mode files (in bytes): %d'%mode_size
print 'Total size of svd files (in bytes): %d'%svd_size
print 'SVD compression ratio: %f'%(float(svd_size) / float(mode_size))