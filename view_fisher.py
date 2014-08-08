#!/usr/bin/env python

"""Plot the results of power spectrum estimation.

:Authors: Shifan Zuo
:Date: 2014-04-08
:email: sfzuo@bao.ac.cn
:usage:
    python view_fisher.py [-h] [--xmin XMIN] [--xmax XMAX] [--ymin YMIN] [--ymax YMAX] [-l FIGLENGTH] [-w FIGWIDTH] [forecast_ps] [obs_ps]
"""

import argparse


def plt_fisher(args):
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    # from mpl_toolkits.axes_grid import AxesGrid

    with h5py.File(args.forecast_ps, 'r') as f1, h5py.File(args.obs_ps, 'r') as f2:
        if f1.attrs['bandtype'] == 'polar':
            polar = True
            K_bands = f1['k_bands'][...]
            k_start = f1['k_start'][...]
            k_end = f1['k_end'][...]
            k_center = f1['k_center'][...]
            theta_bands = f1['theta_bands'][...]
            theta_start = f1['theta_start'][...]
            theta_end = f1['theta_end'][...]
            theta_center = f1['theta_center'][...]
        elif f1.attrs['bandtype'] == 'cartesian':
            polar = False
            kpar_bands = f1['kpar_bands'][...]
            kpar_start = np.unique(np.sort(f1['kpar_start'][...]))
            kpar_end = np.unique(np.sort(f1['kpar_end'][...]))
            kpar_center = np.unique(np.sort(f1['kpar_center'][...]))
            kperp_bands = f1['kperp_bands'][...]
            kperp_start = np.unique(np.sort(f1['kperp_start'][...]))
            kperp_end = np.unique(np.sort(f1['kperp_end'][...]))
            kperp_center = np.unique(np.sort(f1['kperp_center'][...]))
        else:
            raise Exception('Unsupported band type %s.' % f1.attrs['bandtype'])

        forecast_error = f1['errors'][...]
        forecast_fisher = f1['fisher'][...]
        forecast_covariance = f1['covariance'][...]
        forecast_correlation = f1['correlation'][...]
        obs_error = f2['error'][...]
        obs_fisher = f2['fisher'][...]
        obs_covariance = f2['covariance'][...]
        obs_correlation = f2['correlation'][...]
        powerspectrum = f2['powerspectrum'][...]
        
    # plt.figure(figsize=(args.figlength, args.figwidth))
    # fig, axes = plt.subplots(nrows=1, ncols=2)
    # # vmin = min(np.min(np.abs(forecast_fisher)), np.min(np.abs(obs_fisher)))
    # # vmax = max(np.max(np.abs(forecast_fisher)), np.max(np.abs(obs_fisher)))
    # # im = axes.flat[0].imshow(np.log10(np.abs(forecast_fisher)), vmin=vmin, vmax=vmax, origin='lower')
    # # # axes.flat[0].title('Forecasted Fisher matrix')
    # # im = axes.flat[1].imshow(np.log10(np.abs(obs_fisher)), vmin=vmin, vmax=vmax, origin='lower')
    # # # axes.flat[0].title('Observed Fisher matrix')
    # # fig.colorbar(im, shrink=0.75)
    # for data, ax, idx in zip((np.log10(np.abs(forecast_fisher)), np.log10(np.abs(obs_fisher))), axes.flat, (1, 2)):
    #     im = ax.imshow(data, origin='lower')
    #     plt.title('fisher%d' % idx)
    #     # ax.title('fisher')
    # # Make an axis for the colorbar on the right side
    # # plt.title('fisher')
    # cax = fig.add_axes([0.92, 0.2, 0.025, 0.5])
    # fig.colorbar(im, cax=cax)
    # plt.savefig('fisher_matrix.png')

    # fig=plt.figure(figsize=(args.figlength, args.figwidth))
    # grid = AxesGrid(fig, 111, nrows_ncols=(1, 2),
    #             axes_pad = 0.2,
    #             share_all=True,
    #             label_mode = "L",
    #             cbar_location = "right",
    #             cbar_mode="single",
    #             cbar_size='5%'
    #             )
    
    # im = grid[0].imshow(np.log10(np.abs(obs_fisher)))
    # # im.title('fisher')
    # grid.cbar_axes[0].colorbar(im)
    # im = grid[1].imshow(np.log10(np.abs(forecast_fisher)))
    # # im = grid[2].imshow(np.random.random((10,50)))
    # # pl.show()
    # plt.title('fisher')
    # fig.savefig('fisher_matrix.png')

    # fig = plt.figure(figsize=(args.figlength, args.figwidth))
    # plt.subplot(121)
    # im = plt.imshow(np.log10(np.abs(forecast_covariance)), origin='lower')
    # plt.title('Forecasted covariance matrix')
    # # fig.colorbar(shrink=0.75)
    # plt.subplot(122)
    # im = plt.imshow(np.log10(np.abs(obs_covariance)), origin='lower')
    # plt.title('Observed covariance matrix')
    # fig.colorbar(im, shrink=0.75)
    # fig.savefig('covariance_matix.png')

    if polar:
        # Plot the Fisher matrix
        plt.figure(figsize=(args.figlength, args.figwidth))
        plt.subplot(121)
        # plt.imshow(np.log10(np.abs(forecast_fisher)), origin='lower')
        plt.pcolor(np.log10(np.abs(forecast_fisher)))
        plt.title('Forecasted Fisher matrix')
        plt.colorbar()
        plt.subplot(122)
        # plt.imshow(np.log10(np.abs(obs_fisher)), origin='lower')
        plt.pcolor(np.log10(np.abs(obs_fisher)))
        plt.title('Observed Fisher matrix')
        plt.colorbar()
        plt.savefig('fisher_matrix.png')
     
        # Plot covariance matrix
        plt.figure(figsize=(args.figlength, args.figwidth))
        plt.subplot(121)
        # plt.imshow(np.log10(np.abs(forecast_covariance)), origin='lower')
        plt.pcolor(np.log10(np.abs(forecast_covariance)))
        plt.title('Forecasted covariance matrix')
        plt.colorbar()
        plt.subplot(122)
        # plt.imshow(np.log10(np.abs(obs_covariance)), origin='lower')
        plt.pcolor(np.log10(np.abs(obs_covariance)))
        plt.title('Observed covariance matrix')
        plt.colorbar()
        plt.savefig('covariance_matrix.png')
     
        # Plot correlation matrix
        plt.figure(figsize=(args.figlength, args.figwidth))
        plt.subplot(121)
        # plt.imshow(np.log10(np.abs(forecast_correlation)), origin='lower')
        plt.pcolor(np.log10(np.abs(forecast_correlation)))
        plt.title('Forecasted correlation matrix')
        plt.colorbar()
        plt.subplot(122)
        # plt.imshow(np.log10(np.abs(obs_correlation)), origin='lower')
        plt.pcolor(np.log10(np.abs(obs_correlation)))
        plt.title('Observed correlation matrix')
        plt.colorbar()
        plt.savefig('correlation_matrix.png')
     
        # Power spectrum data
        # powerspectrum = np.abs(powerspectrum)
        delta_k = min(k_center[1:] - k_center[:-1]) # spacing between two adjacent k_center
        k_center_offset = k_center + 0.3 * delta_k
        obs_error1 = np.ndarray(shape=(2, len(obs_error)))
        obs_error1[0] = obs_error
        obs_error1[1] = obs_error
        forecast_error1 = np.ndarray(shape=(2, len(forecast_error)))
        forecast_error1[0] = forecast_error
        forecast_error1[1] = forecast_error
        for idx in range(len(obs_error)):
            if obs_error1[0][idx] >= powerspectrum[idx]:
                obs_error1[0][idx] = powerspectrum[idx] - 1.0e-8
            if forecast_error1[0][idx] >= powerspectrum[idx]:
                forecast_error1[0][idx] = powerspectrum[idx] - 1.0e-8
                
        # Plot power spectrum
        # plt.figure(figsize=(args.figlength, args.figwidth))
        plt.figure(figsize=(8, 6))
        plt.errorbar(k_center, powerspectrum, yerr=obs_error1, fmt='bo', ecolor='r', label='Obseved errors')
        plt.errorbar(k_center_offset, powerspectrum, yerr=forecast_error1, fmt='bo', ecolor='g', label='Forcasted errors')
        plt.loglog()
        xmin = args.xmin or np.min(k_center)
        xmax = args.xmax or np.max(k_center)
        plt.xlim(xmin=xmin, xmax=xmax)
        plt.ylim(ymin=args.ymin, ymax=args.ymax)
        plt.xlabel('$k \ [h \ Mpc^{-1}]$')
        plt.ylabel('$P(k) \ [h^{-3}Mpc^3]$')
        plt.legend(numpoints=1)
        plt.savefig('power_spectrum.png')

    else:
        powerspectrum = powerspectrum.reshape(kperp_center.shape[0], -1)
        plt.figure(figsize=(8, 6))
        plt.pcolor(kperp_center, kpar_center, powerspectrum, vmin=-1e2, vmax=1e2)
        plt.xlabel('$k_\perp \ / \ h \ Mpc^{-1}$')
        plt.ylabel('$k_\parallel \ / h \ Mpc^{-1}$')
        plt.colorbar()
        plt.savefig('power_spectrum_2D.png')


parser = argparse.ArgumentParser(description='Plot the results of power spectrum estimation.')
parser.add_argument('forecast_ps', type=str, nargs='?', help='Input hdf5 forecasted fisher file.')
parser.add_argument('obs_ps', type=str, nargs='?', default='', help='Input hdf5 observation fisher file.')
# parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file (png/eps) to save into. If not present, the output image file name will be auto created from the first input map file name and input args (png).')
parser.add_argument('--xmin', type=float, help='Min value of x-axes to plot.')
parser.add_argument('--xmax', type=float, help='Max value of x-axes to plot.')
parser.add_argument('--ymin', type=float, default=1.0e3, help='Min value of y-axes to plot.')
parser.add_argument('--ymax', type=float, default=1.0e5, help='Max value of y-axes to plot.')
parser.add_argument('-l', '--figlength', type=float, default=13, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=5, help='Output figure width.')
# parser.add_argument('-g', '--grid', action='store_true', help='Add meridians and parallels.')
parser.set_defaults(func=plt_fisher)

args = parser.parse_args()
args.func(args)






# import h5py

# import os

# import matplotlib as mpl
# mpl.use('Agg')

# import matplotlib.pyplot as plt

# import numpy

# def psread(psdir):
#    psf = h5py.File(psdir, 'r')
#    fisher_m = psf['fisher_m']
#    fisher_all = psf['fisher_all']
#    #print fisher_all.shape
#    #print fisher_all.value

#    fisher_data = fisher_all.value.real
#    #plt.imshow(numpy.log10(fisher_data), origin='lower' )
#    #plt.colorbar()
#    ##plt.show()
#    #plt.savefig('fisher.png', format='png')


#    c = numpy.linalg.inv(fisher_data)
#    e = numpy.sqrt(c.diagonal())
#    print e

#    #plt.imshow(numpy.log10(c), origin='lower')
#    #plt.colorbar()
#    #plt.savefig('fisher_inv.png', format='png')

#    bandpower = psf['bandpower']
#    bandcenter = psf['bandcenter']
#    p = bandpower.value
#    k = bandcenter.value

#    print p
#    print k

#    error = numpy.ndarray(shape=(2,len(e)))
#    error[0] = e
#    error[1] = e
#    for ii in range(len(e)):
#       if error[0][ii]>p[ii]:
#          error[0][ii]=p[ii]-1.e-10
#    plt.errorbar(k, p, error, fmt='o', c='0.5', capsize=4.5, elinewidth=2)
#    plt.loglog()
#    plt.ylim(ymin=1.e-1, ymax=1.e9)
#    plt.savefig('ps.png', format='png')

#    psf.close()


# if __name__=="__main__":
#    os.environ['SCRATCH'] = '/home/zuoshifan/programming/python/simulation_Richard/myworkspace'
#    scratch = os.getenv("SCRATCH")
#    #print scratch
#    # psread( scratch + "cylinder/voltest/ev/ps/fisher.hdf5" )
#    psread( scratch + "/cylinder/voltest/ev/ps/fisher.hdf5" )
