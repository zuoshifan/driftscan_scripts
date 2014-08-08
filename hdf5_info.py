#!/usr/bin/env python

"""View informations of a hdf5 data file.

:Authors: Shifan Zuo
:Date: 2014-04-02
:email: sfzuo@bao.ac.cn
:usage: python hdf5_info.py [-h] h5files [h5files ...]
"""

import argparse


def h5_info(args):
    """View informations of a hdf5 data file.
    """
    import h5py

    def print_info(name, obj):
        try:
            shape = obj.shape # 'Group' object has no attribute 'shape'
            print name, '  shape = ', shape
        except:
            print name
        # print group/dataset attributes
        for name, value in obj.attrs.iteritems():
	    print name + ":", value

    for h5file in args.h5files:
        print 'File: ', h5file
        with h5py.File(h5file, 'r') as f:
            # print file attributes
            for name, value in f.attrs.iteritems():
                print name + ':', value
            f.visititems(print_info)
        print '-' * 60

parser = argparse.ArgumentParser(description='View informations of a hdf5 data file.')
parser.add_argument('h5files', type=str, nargs='+', help='Input hdf5 files.')
parser.set_defaults(func=h5_info)

args = parser.parse_args()
args.func(args)