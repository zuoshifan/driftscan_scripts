#!/usr/bin/env python

"""view a cartesian map after cut for some `u`, `v` fourier mode.

:Authors: Shifan Zuo
:Date: 2014-08-26
:email: sfzuo@bao.ac.cn
:usage:
    python 
"""

import argparse


def visualize_fft2(args):
    """Visualize sky maps in hdf5 files.

    Arguments
    ---------
    args : argparse namespace.
    """
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.infile, 'r') as f:
        fft2_data = f['fft'][...]
    sffts_data = np.fft.fftshift(fft2_data)
    sffts_data = sffts_data.T # transpose to make u along RA direction, v along  DEC direction
    print 'shape (u, v) = ', sffts_data.shape
    zero_u = sffts_data.shape[0] / 2
    zero_v = sffts_data.shape[1] / 2
    max_u = sffts_data.shape[0] - zero_u - 1
    max_v = sffts_data.shape[1] - zero_v - 1
    cut_data = np.zeros_like(sffts_data)
    umin = max(args.umin, 0)
    umax = min(args.umax, max_u) if args.umax is not None else max_u
    vmin = max(args.vmin, 0)
    vmax = min(args.vmax, max_v) if args.vmax is not None else max_v
    cut_data[(zero_u + umin):(zero_u + umax + 1), (zero_v + vmin):(zero_v + vmax + 1)] = sffts_data[(zero_u + umin):(zero_u + umax + 1), (zero_v + vmin):(zero_v + vmax + 1)]
    cut_data[(zero_u - umax):(zero_u - umin + 1), (zero_v - vmax):(zero_v - vmin + 1)] = sffts_data[(zero_u - umax):(zero_u - umin + 1), (zero_v - vmax):(zero_v - vmin + 1)] # NOTE: in order to have symmetric region, indeces here not symmetric with above
    # special tackle with even case
    if sffts_data.shape[0] % 2 == 0 and args.umax is None:
        cut_data[0, :] = sffts_data[0, :]
    if sffts_data.shape[1] % 2 == 0 and args.vmax is None:
        cut_data[:, 0] = sffts_data[:, 0]
    cut_data = cut_data.T # transpose for we have transposed sffts_data before
    is_cut = np.fft.ifftshift(cut_data)
    cut_map = np.fft.ifft2(is_cut)
    cut_map = cut_map.real   # real sky temperature

    out_file = args.outfile or 'cutuv_%d_%d_%d_%d.' % (umin, umax, vmin, vmax) + args.figfmt
    # plot the cut map
    plt.figure(figsize=(args.figlength, args.figwidth))
    ext = [285, 255, -45, -15]
    plt.imshow(cut_map, extent=ext, origin='lower', vmin=args.min, vmax=args.max)
    plt.xlabel('RA / deg')
    plt.ylabel('DEC / deg')
    plt.colorbar()
    plt.savefig(out_file)



parser = argparse.ArgumentParser(description='view a cartesian map after cut for some `u`, `v` fourier mode.')
parser.add_argument('infile', type=str, nargs='?', help='Input hdf5 fourier transform of a cartensian sky map file.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into. If not present, the output image file name will be auto created from the first input map file name and input args.')
parser.add_argument('-f', '--figfmt', default='pdf', help='Output image format.')
parser.add_argument('--umin', type=int, nargs='?', default=0, help='Min cut threshold of u.')
parser.add_argument('--umax', type=int, nargs='?', help='Max cut threshold of u.')
parser.add_argument('--vmin', type=int, nargs='?', default=0, help='Min cut threshold of v.')
parser.add_argument('--vmax', type=int, nargs='?', help='Max cut threshold of v.')
parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.add_argument('-g', '--grid', action='store_false', help='Add meridians and parallels.')
parser.set_defaults(func=visualize_fft2)

args = parser.parse_args()
args.func(args)