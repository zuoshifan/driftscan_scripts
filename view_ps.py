#!/usr/bin/env python

"""Plot the results of power spectrum estimation.

:Authors: Shifan Zuo
:Date: 2014-04-08
:email: sfzuo@bao.ac.cn
:usage:
    python view_ps.py [-h] [--xmin XMIN] [--xmax XMAX] [--ymin YMIN] [--ymax YMAX] [--vmin VMIN] [--vmax VMAX] [-l FIGLENGTH] [-w FIGWIDTH] [--log] [ps_file]
"""

import argparse
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
# from mpl_toolkits.axes_grid import AxesGrid



def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def matrix_plot(data, figname='fig', figfmt='png', figsize=(8, 6), log=False):
    # mask nan or other invalid values
    data = np.ma.masked_invalid(data, copy=False)
    plt.figure(figsize=figsize)
    if log:
        # plt.pcolor(np.log10(np.abs(data)))
        data = np.ma.masked_where(data <= 0.0, data)
        plt.pcolor(np.log10(data))
    else:
        plt.pcolor(data)
    plt.xlim(xmin=0, xmax=data.shape[1])
    plt.ylim(ymin=0, ymax=data.shape[0])
    plt.colorbar()
    plt.savefig(figname + '.' + figfmt)

def wf_plot(data, figname='wf', figfmt='png', figsize=(13, 5)):
    # mask nan or other invalid values
    data = np.ma.masked_invalid(data, copy=False)
    plt.figure(figsize=(8, 6))
    for row in range(data.shape[0]):
        plt.plot(data[row, :])
    plt.xlim(xmin=0, xmax=data.shape[1])
    plt.xlabel('N')
    plt.ylabel('Window function')
    plt.savefig(figname + '.' + figfmt)

def ps_1d_plot(band_k_center, band_power, k_center, powerspectrum, error, line_type='bo', ps_type='Unwindowed', figname='ps', figfmt='png', figsize=(8, 6), xmin=None, xmax=None, ymin=1.0, ymax=1.0e6):
    if not type(powerspectrum) is list and not type(powerspectrum) is tuple:
        powerspectrum = [powerspectrum]
    if not type(error) is list and not type(error) is tuple:
        error = [error]
    if not type(line_type) is list and not type(line_type) is tuple:
        line_type = [line_type]
    if not type(ps_type) is list and not type(ps_type) is tuple:
        ps_type = [ps_type]
    if len(powerspectrum) != len(error) and len(powerspectrum) != len(ps_type) and len(powerspectrum) != len(line_type):
        raise Exception('Powerspectrum, error, line_type and ps_type do not align.')

    plt.figure(figsize=figsize)
    plt.plot(band_k_center, band_power, 'k-', linewidth=2.5, label = 'Theoretic power spectrum')
    for (ps, err, lt, pt) in zip(powerspectrum, error, line_type, ps_type):
        err1 = np.ndarray(shape=(2, len(err)))
        err1[0] = err
        err1[1] = err
        for idx in range(len(err)):
            if err[idx] >= ps[idx]:
                err1[0][idx] = ps[idx] - 1.0e-8 # err1[0] is lower errorbar
        # plt.errorbar(k_center, ps, yerr=error1, fmt='o', ecolor='r', label='Power spectrum and errors')
        plt.errorbar(k_center, ps, yerr=err1, fmt=lt, label=pt + ' power spectrum')
    plt.loglog()
    xmin = xmin if xmin is not None else np.min(k_center)
    xmax = xmax if xmax is not None else np.max(k_center)
    print 'xmin = %f' % xmin
    print 'xmax = %f' % xmax
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xlabel('$k \ [h \ Mpc^{-1}]$')
    plt.ylabel('$P(k) \ [K^2\ h^{-3}Mpc^3]$')
    plt.legend(numpoints=1)
    plt.savefig(figname + '.' + figfmt)

def ps_2d_plot(kperp_bands, kpar_bands, powerspectrum, figname='ps2d', figfmt='png', figsize=(8, 6), log=False, vmin=None, vmax=None):
    ps2d = powerspectrum.reshape(kperp_bands.shape[0] - 1, kpar_bands.shape[0] - 1).T
    # mask nan or other invalid values
    ps2d = np.ma.masked_invalid(ps2d, copy=False)
    vmin = vmin if vmin is not None else np.min(ps2d)
    vmax = vmax if vmax is not None else np.max(ps2d)
    plt.figure(figsize=figsize)
    if log:
        plt.pcolor(kperp_bands, kpar_bands, np.log10(ps2d), vmin=vmin, vmax=vmax)
    else:
        plt.pcolor(kperp_bands, kpar_bands, ps2d, vmin=vmin, vmax=vmax)
    plt.xlim(xmin=kperp_bands[0], xmax=kperp_bands[-1])
    plt.ylim(ymin=kpar_bands[0], ymax=kpar_bands[-1])
    plt.xlabel('$k_\perp \ / \ h \ Mpc^{-1}$')
    plt.ylabel('$k_\parallel \ / \ h \ Mpc^{-1}$')
    plt.colorbar()
    plt.savefig(figname + '.' + figfmt)



def plt_fisher(args):
    with h5py.File(args.ps_file, 'r') as f:
        if f.attrs['bandtype'] == 'polar':
            polar = True
            K_bands = f['k_bands'][...]
            k_start = f['k_start'][...]
            k_end = f['k_end'][...]
            k_center = f['k_center'][...]
            theta_bands = f['theta_bands'][...]
            theta_start = f['theta_start'][...]
            theta_end = f['theta_end'][...]
            theta_center = f['theta_center'][...]
        elif f.attrs['bandtype'] == 'cartesian':
            polar = False
            kpar_bands = f['kpar_bands'][...]
            kpar_start = np.unique(np.sort(f['kpar_start'][...]))
            kpar_end = np.unique(np.sort(f['kpar_end'][...]))
            kpar_center = np.unique(np.sort(f['kpar_center'][...]))
            kperp_bands = f['kperp_bands'][...]
            kperp_start = np.unique(np.sort(f['kperp_start'][...]))
            kperp_end = np.unique(np.sort(f['kperp_end'][...]))
            kperp_center = np.unique(np.sort(f['kperp_center'][...]))
        else:
            raise Exception('Unsupported band type %s.' % f1.attrs['bandtype'])

        band_power = f['band_power'][...]
        fisher = f['fisher'][...]
        bias = f['bias'][...]

        # Unwindowed
        uw_cv = f['uw_covariance'][...]
        uw_cr = f['uw_correlation'][...]
        uw_wf = f['uw_wfunction'][...]
        uw_err = f['uw_errors'][...]
        uw_ps = f['uw_powerspectrum'][...]

        # Uncorrelated
        uc_cv = f['uc_covariance'][...]
        uc_cr = f['uc_correlation'][...]
        uc_wf = f['uc_wfunction'][...]
        uc_err = f['uc_errors'][...]
        uc_ps = f['uc_powerspectrum'][...]

        # Minimum variance
        mv_cv = f['mv_covariance'][...]
        mv_cr = f['mv_correlation'][...]
        mv_wf = f['mv_wfunction'][...]
        mv_err = f['mv_errors'][...]
        mv_ps = f['mv_powerspectrum'][...]

        # Inverse variance
        iv_cv = f['iv_covariance'][...]
        iv_cr = f['iv_correlation'][...]
        iv_wf = f['iv_wfunction'][...]
        iv_err = f['iv_errors'][...]
        iv_ps = f['iv_powerspectrum'][...]


    # directory and file name
    fisher_name = 'fisher_matrix'
    uw_dir = './unwindowed/'
    uc_dir = './uncorrelated/'
    mv_dir = './minimum_variance/'
    iv_dir = './inverse_variance/'
    cv_name = 'covariance_matrix'
    cr_name = 'correlation_matrix'
    wf_matrix_name = 'window_matrix'
    log_type = '_log'
    lin_type = '_lin'
    wf_curve_name = 'window_curve'
    ps_name = 'powerspectrum'
    ps2d_name = 'powerspectrum2d'
    bandps2d_name = 'band_power'
    err2d_name = 'error2d'
    tps1d_file = os.path.dirname(__file__) + '/21cm_ps1d.hdf5' # Theoretic 21cm power spectrum data

    # create dirs
    mkdir(uw_dir)
    mkdir(uc_dir)
    mkdir(mv_dir)
    mkdir(iv_dir)

    # Plot the Fisher matrix
    matrix_plot(fisher, figname=fisher_name + log_type, log=True)
    matrix_plot(fisher, figname=fisher_name + lin_type, log=False)

    # Plot covariance matrix
    # Unwindowed
    matrix_plot(uw_cv, figname=uw_dir + cv_name + log_type, log=True)
    matrix_plot(uw_cv, figname=uw_dir + cv_name + lin_type, log=False)
    # Uncorrelated
    matrix_plot(uc_cv, figname=uc_dir + cv_name + log_type, log=True)
    matrix_plot(uc_cv, figname=uc_dir + cv_name + lin_type, log=False)
    # Minimum variance
    matrix_plot(mv_cv, figname=mv_dir + cv_name + log_type, log=True)
    matrix_plot(mv_cv, figname=mv_dir + cv_name + lin_type, log=False)
    # Inverse variance
    matrix_plot(iv_cv, figname=iv_dir + cv_name + log_type, log=True)
    matrix_plot(iv_cv, figname=iv_dir + cv_name + lin_type, log=False)

    # Plot correlation matrix
    # Unwindowed
    matrix_plot(uw_cr, figname=uw_dir + cr_name + log_type, log=True)
    matrix_plot(uw_cr, figname=uw_dir + cr_name + lin_type, log=False)
    # Uncorrelated
    matrix_plot(uc_cr, figname=uc_dir + cr_name + log_type, log=True)
    matrix_plot(uc_cr, figname=uc_dir + cr_name + lin_type, log=False)
    # Minimum variance
    matrix_plot(mv_cr, figname=mv_dir + cr_name + log_type, log=True)
    matrix_plot(mv_cr, figname=mv_dir + cr_name + lin_type, log=False)
    # Inverse variance
    matrix_plot(iv_cr, figname=iv_dir + cr_name + log_type, log=True)
    matrix_plot(iv_cr, figname=iv_dir + cr_name + lin_type, log=False)

    # Plot window function matrix
    # Unwindowed
    matrix_plot(uw_wf, figname=uw_dir + wf_matrix_name + log_type, log=True)
    matrix_plot(uw_wf, figname=uw_dir + wf_matrix_name + lin_type, log=False)
    # Uncorrelated
    matrix_plot(uc_wf, figname=uc_dir + wf_matrix_name + log_type, log=True)
    matrix_plot(uc_wf, figname=uc_dir + wf_matrix_name + lin_type, log=False)
    # Minimum variance
    matrix_plot(mv_wf, figname=mv_dir + wf_matrix_name + log_type, log=True)
    matrix_plot(mv_wf, figname=mv_dir + wf_matrix_name + lin_type, log=False)
    # Inverse variance
    matrix_plot(iv_wf, figname=iv_dir + wf_matrix_name + log_type, log=True)
    matrix_plot(iv_wf, figname=iv_dir + wf_matrix_name + lin_type, log=False)

    # Plot window function curve
    # Unwindowed
    # delta function, ignore
    # Uncorrelated
    wf_plot(uc_wf, figname=uc_dir + wf_curve_name)
    # Minimum variance
    wf_plot(mv_wf, figname=mv_dir + wf_curve_name)
    # Inverse variance
    wf_plot(iv_wf, figname=iv_dir + wf_curve_name)


    # Read in theoretic 1d power spectrum data
    with h5py.File(tps1d_file, 'r') as f:
        tk_center = f['k_center'][...]
        tps1d = f['powerspectrum'][...]

    # Plot power spectrum
    if polar:
        # ps = [uw_ps, uc_ps, mv_ps, iv_ps]
        # err = [uw_err, uc_err, mv_err, iv_err]
        # lt = ['bo', 'gs', 'r*', 'c^']
        # pt = ['Unwindowed', 'Uncorrelated', 'Minimum variance', 'Inverse variance']
        ps = [uw_ps, uc_ps, mv_ps]
        err = [uw_err, uc_err, mv_err]
        lt = ['bo', 'gs', 'r*']
        pt = ['Unwindowed', 'Uncorrelated', 'Minimum variance']
        # Unwindowed
        ps_1d_plot(tk_center, tps1d, k_center, ps[0], err[0], line_type=lt[0], ps_type=pt[0], figname=uw_dir + ps_name, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax)
        # Uncorrelated
        ps_1d_plot(tk_center, tps1d, k_center, ps[1], err[1], line_type=lt[1], ps_type=pt[1], figname=uc_dir + ps_name, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax)
        # Minimum variance
        ps_1d_plot(tk_center, tps1d, k_center, ps[2], err[2], line_type=lt[2], ps_type=pt[2], figname=mv_dir + ps_name, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax)
        # Inverse variance
        # ps_1d_plot(tk_center, tps1d, k_center, ps[3], err[3], line_type=lt[3], ps_type=pt[3], figname=iv_dir + ps_name, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax)
        # All together
        ps_1d_plot(tk_center, tps1d, k_center, ps, err, line_type=lt, ps_type=pt, figname=ps_name, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax)

    else:
        # Common band power
        ps_2d_plot(kperp_bands, kpar_bands, band_power, figname=bandps2d_name, vmin=args.vmin, vmax=args.vmax)
        # Unwindowed
        ps_2d_plot(kperp_bands, kpar_bands, uw_ps, figname=uw_dir + ps2d_name, vmin=args.vmin, vmax=args.vmax)
        ps_2d_plot(kperp_bands, kpar_bands, uw_err, figname=uw_dir + err2d_name, vmin=args.vmin, vmax=args.vmax)
        # Uncorrelated
        ps_2d_plot(kperp_bands, kpar_bands, uc_ps, figname=uc_dir + ps2d_name, vmin=args.vmin, vmax=args.vmax)
        ps_2d_plot(kperp_bands, kpar_bands, uc_err, figname=uc_dir + err2d_name, vmin=args.vmin, vmax=args.vmax)
        # Minimum variance
        ps_2d_plot(kperp_bands, kpar_bands, mv_ps, figname=mv_dir + ps2d_name, vmin=args.vmin, vmax=args.vmax)
        ps_2d_plot(kperp_bands, kpar_bands, mv_err, figname=mv_dir + err2d_name, vmin=args.vmin, vmax=args.vmax)
        # Inverse variance
        ps_2d_plot(kperp_bands, kpar_bands, iv_ps, figname=iv_dir + ps2d_name, vmin=args.vmin, vmax=args.vmax)
        ps_2d_plot(kperp_bands, kpar_bands, iv_err, figname=iv_dir + err2d_name, vmin=args.vmin, vmax=args.vmax)




parser = argparse.ArgumentParser(description='Plot the results of power spectrum estimation.')
parser.add_argument('ps_file', type=str, nargs='?', help='Input hdf5 power spectrum file.')
# parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file (png/eps) to save into. If not present, the output image file name will be auto created from the first input map file name and input args (png).')
parser.add_argument('--xmin', type=float, help='Min value of x-axes to plot.')
parser.add_argument('--xmax', type=float, help='Max value of x-axes to plot.')
parser.add_argument('--ymin', type=float, default=1.0, help='Min value of y-axes to plot.')
parser.add_argument('--ymax', type=float, default=1.0e6, help='Max value of y-axes to plot.')
parser.add_argument('--vmin', type=float, default=None, help='Min value of colorbar to plot.')
parser.add_argument('--vmax', type=float, default=None, help='Max value of colorbar to plot.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
# parser.add_argument('-g', '--grid', action='store_true', help='Add meridians and parallels.')
parser.add_argument('--log', action="store_true", help='Use log color scale if present.')
parser.set_defaults(func=plt_fisher)

args = parser.parse_args()
args.func(args)
