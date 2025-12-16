import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils
import os
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as const
import astropy.units as u
import statmorph
from astropy.visualization import simple_norm
from statmorph.utils.image_diagnostics import make_figure
from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import detect_sources
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.segmentation import deblend_sources
from photutils.segmentation import SourceFinder
from photutils.segmentation import SourceCatalog
from astropy.visualization import LogStretch
from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats, CircularAperture
from photutils.datasets import make_4gaussians_image
from photutils.aperture import aperture_photometry
from photutils.isophote import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import simple_norm
from scipy.optimize import leastsq
from scipy.optimize import fsolve
from scipy.optimize import minimize
import matplotlib.ticker as ticker
import copy
import sys
import os
import math
import logging
# import galsim

plt.rcParams['figure.figsize'] = [16, 20]
plt.rcParams.update({'font.size': 15})
plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.interpolation'] = 'nearest'
plt.rcParams['image.cmap'] = 'viridis'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.style"] = "normal"

def fetch_filename(band, names, idx, sept):
    """ 
    fetch the filename from the list of filenames with given criteria 
    where
        num is the number index of the file to be fetched
        names is the list of filenames
        idx is the index of the str object 'num' in the filename 
        sept is the separator of the filename
    """
    for i, t in enumerate(names):
        if band == t.split(sept)[idx].upper():
            return t


def hsc_filter(band, path):
    """
    return the throughput of HSC filters derived from the SOV website
    """
    return np.genfromtxt(path + 'Subaru_HSC.{}_filter.dat'.format(band))


def abmag2flambda(mag, lambda_m):
    """
    convert AB magnitude to f_lambda
    where
        mag is the AB magnitude
        lambda_m is the wavelength in micron
    """
    f_nu = 10**(-0.4*(mag+48.6)) * u.Unit("erg s-1 cm-2 Hz-1")
    lambda_m = lambda_m * u.AA
    f_lambda = (const.c / lambda_m**2) * f_nu
    return f_lambda.to("erg s-1 cm-2 AA-1")


def rad2deg(rad):
    """
    convert radians to degrees
    note that for galfitm, the position angle is measured from the positive y-axis where
        position angle (PA) [deg: Up=0, Left=90]
        while for the statmorph, the PA is measured from the positive x-axis where
        position angle (PA) [deg: Up=90, Right=0]
    """
    return rad * 180 / np.pi

def deg2rad(deg):
    """
    convert degrees to radians
    note that for photutils, the position angle is measured from the positive x-axis where
        position angle (PA) [deg: Up=90, Right=0]
    """
    return deg * np.pi / 180

def root_mean_sq(arr):
    """
    calculate the root mean square of the array
    """
    return np.sqrt(np.mean(arr**2))

def root_sum_sq(arr):
    """
    calculate the root of the sum of square of the array
    """
    return np.sqrt(np.sum(arr**2))

def elip2ar(elip):
    """
    convert ellipticity to axis ratio (ar)
    """
    return np.sqrt(1 - elip**2)

def ar2elip(ar):
    """
    convert axis ratio (ar) to ellipticity
    """
    return np.sqrt(1 - ar**2)

def area_ellipse(a, ar):
    """
    calculate the area of the ellipse
    where
        a is the semi-major axis
        ar is the axis ratio
    """
    return np.pi * a**2 * ar

def hsc_waves():
    hsc_g_wave = 4669.70        # lambda_effective, and filter only
    hsc_r_wave = 6135.56
    hsc_i_wave = 7656.73
    hsc_z_wave = 8919.85
    hsc_y_wave = 9993.28

    hsc_g_lam_min = 3942.85	
    hsc_g_lam_max = 5547.08
    hsc_r_lam_min = 5330.23	
    hsc_r_lam_max = 7071.80
    hsc_i_lam_min = 6800.59
    hsc_i_lam_max = 8658.53	
    hsc_z_lam_min = 8364.41
    hsc_z_lam_max = 9498.95
    hsc_y_lam_min = 9000.00	
    hsc_y_lam_max = 10934.85
    return np.array([[hsc_g_wave, hsc_r_wave, hsc_i_wave, hsc_z_wave, hsc_y_wave], \
        [hsc_g_lam_min, hsc_r_lam_min, hsc_i_lam_min, hsc_z_lam_min, hsc_y_lam_min], \
        [hsc_g_lam_max, hsc_r_lam_max, hsc_i_lam_max, hsc_z_lam_max, hsc_y_lam_max]])
            


def morph_est(img, psf0, sigma, plot=False):
    """
    estimate the morphology of the galaxy
    where
        img is the image of the galaxy
        psf is the psf of the image
        sigma is the noise level of the image (weightmap)
        mask is the maskmap of the image
    """
    bkg_estimator = MedianBackground()
    kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0

    bkg = Background2D(img, (15,15), filter_size=(5, 5), bkg_estimator=bkg_estimator)
    data = img - bkg.background  # subtract the background

    threshold = 0.8 * bkg.background_rms
    convolved_data = convolve(data, kernel)

    finder = SourceFinder(npixels=4, progress_bar=False)
    segment_lab = finder(convolved_data, threshold)
    labelarr = segment_lab.data
    segmap = copy.deepcopy(labelarr)
    cent_label = labelarr[labelarr.shape[0]//2, labelarr.shape[1]//2]
    if cent_label != 0:
        segment_lab.remove_labels(cent_label)
        segmap[segmap != cent_label] = 0
        segmap[segmap == cent_label] = 1
    mask = segment_lab.data
    mask[mask>0] = 1

    mask_map = np.array(mask).astype(bool)
    source_morphs = statmorph.source_morphology(data, segmap, mask=mask_map, weightmap=sigma, psf=psf0)
    morph = source_morphs[0]

    if plot:
        try:
            fig = make_figure(morph)
            plt.close(fig)
        except:
            print('failed to plot the morphology with morph_est')

    return morph


def gen_galfitm_feedme(feedme_name, img_idx, qso_names, n_gal, bandnames, Band_names, img_path,
        sig_path, psf_path, mask_path, maskfile, result_name, pix_scale, hsc_waves, fit_region = None,
        iter_n=1, pre_res=None, numid_psf=1, num_comp=2):
    """
    generate the feedme file for galfitm
    where
        feedme_name is the path and name of the feedme file
        bandnames is the list of band names
        img_path is the path of the image_data
        sig_path is the path of the sigma_map
        psf_path is the path of the psf
        mask_path is the path of the mask-files
        n_gal is the number of galaxies in the image
        img_idx is the index of the qso image data in the whole qso images
        qso_names is the list of the qso image data names
        result_name is the path+name of the result file
        Band_names is the list of the band names in bold
        pix_scale is the pixel scale of the image
        hsc_waves is the list of the effective wavelengths of the HSC filters
    """
    n_band = len(bandnames)
    with open(feedme_name, 'w+') as f:
        f.write('A) ')
        for s in range(n_band):
            band = bandnames[s]
            fname = fetch_filename(band, qso_names, idx=1, sept='_')
            f.write(img_path + fname)
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('A1) ')
        for s in range(n_band):
            f.write(bandnames[s])
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('A2) ')
        for s in range(n_band):
            f.write(str(hsc_waves[s]))
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')

        f.write('B) ' + result_name + '\n')
        f.write('C) ')
        for s in range(n_band):
            band = bandnames[s]
            fname = fetch_filename(band, qso_names, idx=1, sept='_')
            f.write(sig_path + fname.replace('data', 'err_scale'))
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('D) ')
        for s in range(n_band):
            band = bandnames[s]
            print(s, band)
            # psf_f = str(img_idx + s*n_gal) + '-psf-' + Band_names[s] + '.fits'
            # fname = fetch_filename(band, qso_names, idx=1, sept='_')\
            psf_file = band + '_psf.fits'
            f.write(psf_path +  psf_file)
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('E) 1'+ '\n')
        f.write('F) ')
        for s in range(n_band):
            # fname = fetch_filename(img_idx + s*n_gal, qso_names, idx=0, sept='-')
            fname = maskfile
            #'early_mask.fits'
            f.write(mask_path +  fname)
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('G) ../galfitm_SS.CONSTRAINTS'+ '\n')
        f.write('H) ')
        img0 = fits.getdata(img_path + fname)
        if(fit_region == None):
            f.write('1    '+str(img0.shape[1])+'   1    '+str(img0.shape[0]) + '\n')    # when the 'error and data must have the same shape' error occurs, change the order here.
        else:
            f.write(str(fit_region[0][0]) + '    '+str(fit_region[0][1]) +   '    '+ str(fit_region[1][0]) + '    '+ str(fit_region[1][1]) + '\n')    # when the 'error and data must have the same shape' error occurs, change the order here.
        f.write('I) '+str(img0.shape[1]-4)+'    '+str(img0.shape[0]-4) + ' \n')
        f.write('J) 28.08651939228399,28.08651939228399,28.08651939228399,28.08651939228399,28.08651939228399,28.08651939228399,28.08651939228399,28.08651939228399' + '\n')
        f.write('K) '+str(pix_scale)+'  '+str(pix_scale)+' \n')
        f.write('O) regular' + '\n')
        f.write('P) 0' + '\n')
        f.write('W) input,model,residual,component' + '\n')

        # estimate the initial parameters for fitting input of GalfitM
        # the initial parameters are estimated by the statmorph and photutils
        if iter_n != 1:
            fit_res = fits.open(pre_res)
        num_id_sersic = numid_psf + 1
        if iter_n == 1:
            mags = []
            bkg_estimator = MedianBackground()
            kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0
            for t in range(n_band):
                band = bandnames[s]
                fname = fetch_filename(band, qso_names, idx=1, sept='_')
                img0 = fits.getdata(img_path + fname)
                mask = fits.getdata(mask_path + maskfile)
                img0 *= (1-mask)                # load the mask with central region masked
                bkg = Background2D(img0, (15,15), filter_size=(5, 5), bkg_estimator=bkg_estimator)
                data = img0 - bkg.background    # subtract the background
                threshold = 1.5 * bkg.background_rms
                convolved_data = convolve(data, kernel)
                segment_map = detect_sources(convolved_data, threshold, npixels=10)
                segm_deblend = deblend_sources(convolved_data, segment_map, npixels=10, nlevels=32, contrast=0.001, progress_bar=False)
                cat = SourceCatalog(data, segm_deblend, convolved_data=convolved_data)
                seg_flux = cat.to_table(columns = ['segment_flux'])['segment_flux'][0]
                mag0 = -2.5*np.log10(seg_flux) + 28.08651939228399
                mags.append(mag0)
        else:
            mags = []
            for t in range(n_band):
                mag0_str = fit_res[n_band+t].header[str(numid_psf) + '_MAG_'+Band_names[t]].split(' ')[0]
                mag1_str = fit_res[n_band+t].header[str(num_id_sersic) + '_MAG_'+Band_names[t]].split(' ')[0]
                if ('*' in mag0_str) or ('*' in mag1_str):
                    bkg_estimator = MedianBackground()
                    kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0
                    fname = fetch_filename(img_idx + t*n_gal, qso_names, idx=0, sept='-')
                    img0 = fits.getdata(img_path + fname)
                    mask = fits.getdata(mask_path + 'maskcent-' + fname)
                    img0 *= (1-mask)                # load the mask with central region masked
                    bkg = Background2D(img0, (15,15), filter_size=(5, 5), bkg_estimator=bkg_estimator)
                    data = img0 - bkg.background    # subtract the background
                    threshold = 1.5 * bkg.background_rms
                    convolved_data = convolve(data, kernel)
                    segment_map = detect_sources(convolved_data, threshold, npixels=10)
                    segm_deblend = deblend_sources(convolved_data, segment_map, npixels=10, nlevels=32, contrast=0.001, progress_bar=False)
                    cat = SourceCatalog(data, segm_deblend, convolved_data=convolved_data)
                    seg_flux = cat.to_table(columns = ['segment_flux'])['segment_flux'][0]
                    mag0 = -2.5*np.log10(seg_flux) + 28.08651939228399
                    mag1 = mag0 + 0.1
                    if '*' not in mag0_str:
                        mag0 = float(mag0_str)
                    if '*' not in mag1_str:
                        mag1 = float(mag1_str)
                else:
                    mag0 = float(mag0_str)
                    mag1 = float(mag1_str)
                mags.append(min(mag0, mag1))

        f.write('\n')
        f.write('#estimated mag: ')
        for s in range(n_band):
            f.write(str(mags[s]))
            if s==n_band-1:
                f.write('\n')
            else:
                f.write(',')
        f.write('\n')
        f.write('# -----------------------------------------------------------------------------' + '\n')
        f.write('#   par)    par value(s)    fit toggle(s)    # parameter description' + '\n')
        f.write('# -----------------------------------------------------------------------------' + '\n')

        f.write('\n')

        f.write('# Object number: 1' + '\n')
        f.write(' 0) psf              #  object type' + '\n')
        f.write(' 1) ' + str(img0.shape[1] /2) + '      2       # position x [pixel]' + '\n')
        f.write(' 2) ' + str(img0.shape[0] /2) + '      2       # position y [pixel]' + '\n')
        f.write(' 3) ')
        for s in range(n_band):
            f.write(str(mags[s]))
            if s==n_band-1:
                f.write('      ' + str(n_band) + '       # total magnitude' + '\n')
            else:
                f.write(',')
        f.write(' Z) 0                #  output option (0 = resid., 1 = Don\'t subtract)' + '\n') 

        # estimate morphology of the galaxy with statmorph
        PA = []
        ar = []
        reff = []
        sersic_n = []
        for t in range(n_band):
            band = bandnames[t]
            fname = fetch_filename(band, qso_names, idx=1, sept='_')
            if iter_n == 1:
                img0 = fits.getdata(img_path + fname)
            else:
                numid_psf_model = numid_psf + 3*n_band -1 + t*num_comp
                img0 = fit_res[t].data - fit_res[numid_psf_model].data      # data - psf component
            sig0 = fits.getdata(sig_path + fname.replace('data', 'err'))
            psf_file = band + '_psf.fits'
            psf0 = fits.getdata(psf_path + psf_file)
            # print(img0)
 #            print(psf0)
 #            print(sig0)
            # print(psf_file)
            morph = morph_est(img0, psf0, sig0, plot=True)

            PA.append(rad2deg(morph.orientation_centroid)-90)
            ar.append(elip2ar(morph.sersic_ellip))
            reff.append(morph.sersic_rhalf)
            sersic_n.append(morph.sersic_n)

        sersic_n_mean = np.mean(np.array(sersic_n))
        ar_mean = np.mean(np.array(ar))
        reff_mean = np.mean(np.array(reff))
        PA_mean = np.mean(np.array(PA))
        if sersic_n_mean <= 0 or np.isnan(ar_mean) or reff_mean <= 0:
            if iter_n == 1:
                print(result_name.split('/')[-1] +'  morphology estimation failed for the galaxy for at least one band, something really wrong with the image')
            else:
                print(result_name.split('/')[-1] +'  GalfitM failed to subtract the psf successfully for at least one band')

            if np.array(sersic_n).max() > 0 and  np.sum(np.isnan(ar)) != len(ar)  and np.array(reff).max() > 0:
                sersic_n_mean = np.array(sersic_n).max()
                reff_mean = np.array(reff).max()
                ar_mean = np.array(ar)[~np.isnan(ar)].mean()
            elif iter_n == 1:
                sersic_n_mean = 4
                reff_mean = 5
                ar_mean = 0.8
            else:
                sersic_n_mean = float( fit_res[n_band].header[str(num_id_sersic) + '_N_'+Band_names[0]].split(' ')[0] )
                reff_mean = float( fit_res[n_band].header[str(num_id_sersic) + '_RE_'+Band_names[0]].split(' ')[0] )
                ar_mean = float( fit_res[n_band].header[str(num_id_sersic) + '_AR_'+Band_names[0]].split(' ')[0] )
                print('using the previous result for the morphology estimation from GalfitM')

        f.write('# Object number: 2' + '\n')
        f.write(' 0) sersic             #  object type' + '\n')
        f.write(' 1) ' + str(img0.shape[1] /2) + '      2       # position x [pixel]' + '\n')
        f.write(' 2) ' + str(img0.shape[0] /2) + '      2       # position y [pixel]' + '\n')
        f.write(' 3) ')
        for s in range(n_band):
            f.write(str(mags[s] + 0.3))
            if s==n_band-1:
                f.write('      ' + str(n_band) + '       # total magnitude' + '\n')
            else:
                f.write(',')
        f.write(' 4) '+str(reff_mean)+'      3          #  R_e (half-light radius)   [pix]' + '\n')
        if sersic_n_mean > 8:
            sersic_n_mean = 8
        f.write(' 5) '+str(sersic_n_mean)+'      3          #  Sersic index n (de Vaucouleurs n=4)' + '\n')
        #f.write(' 4) '+str(reff[1])+'      4          #  R_e (half-light radius)   [pix]' + '\n')
        #f.write(' 5) '+str(sersic_n[1])+'      4          #  Sersic index n (de Vaucouleurs n=4)' + '\n')
        f.write(' 6) 0.0000      0          #     -----' + '\n')
        f.write(' 7) 0.0000      0          #     -----' + '\n')
        f.write(' 8) 0.0000      0          #     -----' + '\n')
        f.write(' 9) '+str(ar_mean)+'      3          #  axis ratio (b/a)' + '\n')
        f.write('10) '+str(PA_mean)+'    2          #  position angle (PA) [deg: Up=0, Left=90]' + '\n')
        f.write(' Z) 0                      #  output option (0 = resid., 1 = Don\'t subtract)' + '\n')


def gen_exec_script(script_name, obj_num, feedme_name):
    """
    generate the execution script for galfitm
    where
        script_name is the path+name of the execution script
        obj_num is the number of iterations
        feedme_name is the name of the feedme file
    """
    with open(script_name, 'w+') as f:
        for s in range(obj_num):
            f.write('galfitm ' + feedme_name[s] +'\n')


def res_alys(result_name, numid_psf, num_comp, band_names, Band_names, 
    img_idx, n_gal, qso_names, sig_path, bandwaves,
    apert, pix_scale):
    '''
    numid_psf is the number index of psf in the model
    apert is the aperture radius in arcsec
    '''
    fit_res = fits.open(result_name)
    n_band = len(band_names)
    hd = fit_res[n_band+1].header

    aperture_radius = apert / pix_scale   # arcsec to pixel
    qso_f_lambda = []
    qso_f_lambda_err_up = []
    qso_f_lambda_err_low = []
    gal_f_lambda = []
    gal_f_lambda_err_up = []
    gal_f_lambda_err_low = []
    for t in range(n_band):
        positions = [float(hd[str(numid_psf) + '_XC_'+Band_names[t]].split(' ')[0]), 
            float(hd[str(numid_psf) + '_YC_'+Band_names[t]].split(' ')[0])]
        aperture = CircularAperture(positions, r = aperture_radius)
        fname = fetch_filename(img_idx + t*n_gal, qso_names, idx=0, sept='-')
        sigma_map = fits.getdata(sig_path + 'sig-' + fname)

        numid_psf_model = numid_psf + 3*n_band -1 + t*num_comp
        phot_table = aperture_photometry(fit_res[numid_psf_model].data, 
            aperture, error=sigma_map , method='exact')        # psf component
        flux = phot_table['aperture_sum'][0] 
        mag = -2.5*np.log10(flux) + 28.08651939228399
        f_lam = abmag2flambda(mag, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA

        #flux_err = phot_table['aperture_sum_err'][0]
        #mag_err = -2.5*np.log10(flux_err) + 28.08651939228399
        mag_err = float(hd[str(numid_psf) + '_MAG_'+Band_names[t]].split(' ')[-1])     # mag_err taken as from the galfitm fitting result
        f_lam_err_up = abmag2flambda(mag-mag_err, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA
        f_lam_err_low = abmag2flambda(mag+mag_err, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA

        qso_f_lambda.append(f_lam)
        qso_f_lambda_err_up.append(f_lam_err_up)
        qso_f_lambda_err_low.append(f_lam_err_low)
        #print('flux density in 10^-17 erg/s/cm^2/AA: ', f_lam, f_lam_err)  # show output

    numid_sersic = numid_psf + 1
    for t in range(n_band):
        positions = [float(hd[str(numid_sersic) + '_XC_'+Band_names[t]].split(' ')[0]),
            float(hd[str(numid_sersic) + '_YC_'+Band_names[t]].split(' ')[0])]
        aperture = CircularAperture(positions, r = aperture_radius)
        fname = fetch_filename(img_idx + t*n_gal, qso_names, idx=0, sept='-')
        sigma_map = fits.getdata(sig_path + 'sig-' + fname)

        numid_psf_model = numid_psf + 3*n_band -1 + t*num_comp
        phot_table = aperture_photometry(fit_res[t].data - fit_res[numid_psf_model].data, 
            aperture, error=sigma_map , method='exact')        # data - psf component
        flux = phot_table['aperture_sum'][0] 
        mag = -2.5*np.log10(flux) + 28.08651939228399
        f_lam = abmag2flambda(mag, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA

        #flux_err = phot_table['aperture_sum_err'][0]
        #mag_err = -2.5*np.log10(flux_err) + 28.08651939228399
        mag_err = root_sum_sq(np.array([ float(hd[str(numid_psf) + '_MAG_'+Band_names[t]].split(' ')[-1]),
            float(hd[str(numid_sersic) + '_MAG_'+Band_names[t]].split(' ')[-1]) ]))
            # mag_err for host_gal taken as sum of both the psf and sersic component
        mag_err_qso = float(hd[str(numid_psf) + '_MAG_'+Band_names[t]].split(' ')[-1])
        if mag_err > 100* mag_err_qso:
            mag_err = np.sqrt(2) * mag_err_qso
        f_lam_err_up = abmag2flambda(mag-mag_err, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA
        f_lam_err_low = abmag2flambda(mag+mag_err, bandwaves[t]).value *10**17      # flux density in 10^-17 erg/s/cm^2/AA

        gal_f_lambda.append(f_lam)
        gal_f_lambda_err_up.append(f_lam_err_up)
        gal_f_lambda_err_low.append(f_lam_err_low)
    return np.array([qso_f_lambda, qso_f_lambda_err_up, qso_f_lambda_err_low]), np.array([gal_f_lambda, gal_f_lambda_err_up, gal_f_lambda_err_low])


def plot_galfitm_res(result_name, figure_name, numid_psf, num_comp, band_names,
                     img_idx, img_path, mask_file, n_gal, pixscale, save=False, show=True):
    n_band = len(band_names)
    fig, ax = plt.subplots(n_band,5)        # 5 columns
    fig.set_size_inches(20,18)
    norm = ImageNormalize(stretch=LogStretch(a=10000), vmin=0, vmax=50)
    fig.tight_layout()
    plt.subplots_adjust(wspace = None, hspace=0.1)

    ax[0,0].set_title('data')
    ax[0,1].set_title('model')
    ax[0,2].set_title('data - psf')
    ax[0,3].set_title('residual')
    ax[0,4].set_title('Surface Brightness Profile')

    for t in range(n_band):
        ax[t,0].set_ylabel(band_names[t])

    qso_names = os.listdir(img_path)
    # fname = fetch_filename(img_idx, qso_names, idx=0, sept='-')
    fname = fetch_filename(img_idx, qso_names, idx=1, sept='_')
    img0 = fits.getdata(img_path + fname)
    datasize = img0.shape
    x_cent = (datasize[1]-1)/2
    y_cent = (datasize[0]-1)/2
    for j in range(n_band):
        for k in range(4):
            ax[j,k].set_xticks([x_cent-4/pixscale, x_cent-2/pixscale, x_cent, x_cent+2/pixscale, x_cent+4/pixscale],[-4,-2,0,2,4])
            ax[j,k].set_yticks([y_cent-4/pixscale, y_cent-2/pixscale, y_cent, y_cent+2/pixscale, y_cent+4/pixscale],[-4,-2,0,2,4])
    for j in range(5):
        ax[n_band -1 ,j].set_xlabel('arcsec')

    gal_res = fits.open(result_name)
    hd = gal_res[5].header

    for s in range(n_band):
        im1 = ax[s,0].imshow(gal_res[s].data, norm=norm)        # data
        divider = make_axes_locatable(ax[s,0])
        cax = divider.append_axes('right', size='2%', pad=0.02)
        fig.colorbar(im1, cax=cax, orientation='vertical',format=ticker.LogFormatter())

        im2 = ax[s,1].imshow(gal_res[s + n_band].data, norm=norm)       # model
        divider = make_axes_locatable(ax[s,1])
        cax = divider.append_axes('right', size='2%', pad=0.02)
        fig.colorbar(im2, cax=cax, orientation='vertical',format=ticker.LogFormatter())

        # fname = fetch_filename(img_idx + s*n_gal, qso_names, idx=0, sept='-')   
        band = band_names[s]
        fname = fetch_filename(img_idx, qso_names, idx=1, sept='_')
        
        mask = fits.getdata(mask_file)
        unmask_cent = copy.deepcopy(mask)
        x_cent = mask.shape[0]//2
        y_cent = mask.shape[1]//2
        unmask_cent[x_cent-3:x_cent+3, y_cent-3:y_cent+3] = 0
        numid_psf_model = numid_psf + 3*n_band -1 + s*num_comp
        im3 = ax[s,2].imshow(gal_res[s].data *(1-unmask_cent) - gal_res[numid_psf_model].data, norm=norm)     # masked_data - psf
        divider = make_axes_locatable(ax[s,2])
        cax = divider.append_axes('right', size='2%', pad=0.02)
        fig.colorbar(im3, cax=cax, orientation='vertical',format=ticker.LogFormatter())

        residual = gal_res[s + 2*n_band].data * (1-unmask_cent)
        max_val = max(np.max(residual), -np.min(residual))
        im4 = ax[s,3].imshow(residual, cmap='seismic', origin='lower', vmin=-max_val, vmax=max_val)     # masked_residual
        divider = make_axes_locatable(ax[s,3])
        cax = divider.append_axes('right', size='2%', pad=0.02)
        fig.colorbar(im4, cax=cax, orientation='vertical')

        # numid_sersic = numid_psf + 1
        numid_sersic = 2
        hd = gal_res[s + n_band].header
        ar = float(hd[str(numid_sersic)+ '_AR_'+ str.upper(band_names[s])].split(' ')[0])
        PA = float(hd[str(numid_sersic)+ '_PA_'+ str.upper(band_names[s])].split(' ')[0])
        try:
            reff = float(hd[str(numid_sersic)+ '_RE_'+ str.upper(band_names[s])].split(' ')[0])
        except:
            reff = float(hd[str(numid_sersic)+ '_RE_0'].split(' ')[0])
        positions = [float(hd[str(numid_psf) + '_XC_'+str.upper(band_names[s])].split(' ')[0]), 
            float(hd[str(numid_psf) + '_YC_'+str.upper(band_names[s])].split(' ')[0])]
        geometry = EllipseGeometry(x0=positions[0], y0=positions[1], sma=max(reff, 5.), eps=min(0.4, 1-ar), pa= deg2rad(PA + 90))     
                    # eps = 1 - b/a for photutils, which is actually flattening. PA is (in radians) of the semimajor 
                    # sma ('a') should be at least 2.
                    # axis in relation to the positive x axis of the image array (rotating towards the positive y axis), within the range (0, pi].
        ellipse_data = Ellipse(gal_res[s].data*(1-unmask_cent), geometry)
        isolist_data = ellipse_data.fit_image()
        sb_data = 28.08651939228399- 2.5*np.log10(isolist_data.intens / pixscale**2)
        sb_data_err_low = 28.08651939228399 - 2.5*np.log10((isolist_data.intens + isolist_data.int_err) / pixscale**2)
        sb_data_err_up = 28.08651939228399 - 2.5*np.log10((isolist_data.intens - isolist_data.int_err) / pixscale**2)
        sb_data_err = [sb_data-sb_data_err_low, sb_data_err_up-sb_data]
        ellipse_model = Ellipse(gal_res[s + n_band].data, geometry)
        isolist_model = ellipse_model.fit_image()
        sb_model = 28.08651939228399- 2.5*np.log10(isolist_model.intens / pixscale**2)
        sb_model_err_low = 28.08651939228399 - 2.5*np.log10((isolist_model.intens + isolist_model.int_err) / pixscale**2)
        sb_model_err_up = 28.08651939228399 - 2.5*np.log10((isolist_model.intens - isolist_model.int_err) / pixscale**2)
        sb_model_err = [sb_model-sb_model_err_low, sb_model_err_up-sb_model]
               
        ax[s,4].errorbar(isolist_data.sma * pixscale, sb_data, yerr=sb_data_err, fmt='o', color='red', label='data')
        ax[s,4].plot(isolist_model.sma * pixscale, sb_model, color='blue', label='model')
        ax[s,4].set_xlim(0, 15*pixscale)
        try:
            index = (np.abs(isolist_data.sma - 15)).argmin()
            y_up_cut = 28.08651939228399- 2.5*np.log10(isolist_data.intens[index] / pixscale**2)
            ax[s,4].set_ylim( top = y_up_cut+0.5 )
        except:
            print('failed ellipse fitting of photutils') 
        ax[s,4].invert_yaxis()
        ax[s,4].legend()
    if save:
        fig.savefig(figure_name, dpi=300, bbox_inches='tight')

    plt.close(fig)
    if show:
        plt.show()
