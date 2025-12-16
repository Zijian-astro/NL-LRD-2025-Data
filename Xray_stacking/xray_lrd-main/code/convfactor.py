
# write it for a single object

import numpy as np
from astropy.table import Table

import sherpa.astro.xspec as xspec
# from sherpa.ui import create_model_component
from sherpa.astro.ui import *
def flux_to_luminosity(redshift, flux, nH, mwnH):

    #redshift = 5

    # the model is like this

    create_model_component("powlaw1d", "powlaw")
    powlaw.gamma = 1.8
    powlaw.ref = 1.0
    powlaw.ampl = 1.0

    create_model_component("xsphabs", "absorb_mw")
    create_model_component("xszphabs", "absorb_agn")
    # create_model_component("phabs", "absorb_mw")
    # create_model_component("phabs", "absorb_agn")

    absorb_mw.nH = mwnH

    absorb_agn.redshift=5
    absorb_agn.nH=nH

    srcmdl = powlaw * absorb_mw * absorb_agn

    # what is the flux?
    create_model_component("xscflux", "cflux")
    fluxresult = cflux(srcmdl).calc([3])
    #print(fluxresult.thawedpars[0])

    # we know the trueflux. let's normalize it
    trueflux = 1e-18
    lgtrueflux = np.log10(trueflux)
    factor = 10**(lgtrueflux-fluxresult.thawedpars[0])
    print(factor)

    # define a new model
    create_model_component("powlaw1d", "powlaw2")
    powlaw2.gamma = 2.0
    powlaw2.ref = 1.0
    powlaw2.ampl = factor

    srcmdl2 = powlaw2
    create_model_component("xsclumin", "clumin")

    lum_result = clumin(srcmdl2).calc([3])
    print(lum_result)
    return lum_result.thawedpars[0]

def luminosity_to_flux_factor(redshift, logNH, mwnH):
    create_model_component("powlaw1d", "powlaw")
    powlaw.gamma = 1.8
    powlaw.ref = 1.0
    powlaw.ampl = 100

    create_model_component("xsphabs", "absorb_mw")
    create_model_component("xszphabs", "absorb_agn")
    # create_model_component("phabs", "absorb_mw")
    # create_model_component("phabs", "absorb_agn")
    
    absorb_mw.nH = mwnH/1e22

    absorb_agn.redshift=redshift
    absorb_agn.nH=10**(logNH-22)

    srcmdl = powlaw * absorb_mw * absorb_agn
    # what is the flux?
    create_model_component("xsclumin", "clumin")

    loglum = 42

    clumin.Emin = 2
    clumin.Emax = 10
    clumin.lg10Lum = loglum
    clumin.Redshift = redshift

#    freeze(powlaw.ampl)
#    freeze(powlaw.gamma)

    convsrcmdl = absorb_mw * absorb_agn * clumin(powlaw)

    convsrcmdl_unabs = clumin(powlaw)
    print(convsrcmdl_unabs)

#    xlist = np.arange(0.1, 100, 0.01)
#    ylist = convsrcmdl_unabs(xlist)

    data = dataspace1d(0.1,100,0.01)
    get_data()
    lum_test = 10**loglum

#    set_source(convsrcmdl_unabs)
#    set_analysis('energy')
#    flux_test = calc_energy_flux(0.5, 2, model=convsrcmdl_unabs)
    set_source(convsrcmdl)
    flux_test_soft = calc_energy_flux(0.5, 2, model=convsrcmdl)
    flux_test_hard = calc_energy_flux(2, 8, model=convsrcmdl)
    flux_test_all = calc_energy_flux(0.5, 8, model=convsrcmdl)

    rate_test_soft = calc_photon_flux(0.5, 2, model=convsrcmdl)
    rate_test_hard = calc_photon_flux(2, 8, model=convsrcmdl)
    rate_test_all = calc_photon_flux(0.5, 8, model=convsrcmdl)

    fluxfactor_soft = flux_test_soft / lum_test
    fluxfactor_hard = flux_test_hard / lum_test
    fluxfactor_all = flux_test_all / lum_test

    ratefactor_soft = rate_test_soft / lum_test
    ratefactor_hard = rate_test_hard / lum_test
    ratefactor_all = rate_test_all / lum_test

    return fluxfactor_soft, fluxfactor_hard, fluxfactor_all,\
            ratefactor_soft, ratefactor_hard, ratefactor_all

def get_convert_factor_sample(inname, outname):
    tbl = Table.read(inname)

    for nH in (21,22,23):
        fluxfactor_s = []
        fluxfactor_h = []
        fluxfactor_a = []

        ratefactor_s = []
        ratefactor_h = []
        ratefactor_a = []

        print(nH)
        for index in range(len(tbl)):
            redshift = tbl['redshift'][index]
            mwnH = tbl['NH_MW'][index]
            fs, fh, fA, rs, rh, rA = luminosity_to_flux_factor(redshift,nH, mwnH)
            fluxfactor_s.append(fs)
            fluxfactor_h.append(fh)
            fluxfactor_a.append(fA)

            ratefactor_s.append(rs)
            ratefactor_h.append(rh)
            ratefactor_a.append(rA)

        if np.abs(nH-int(nH))>0.1:
            print('what?')
            nH = 10*nH
        tbl['ff_soft_%d'%nH] = np.array(fluxfactor_s)
        tbl['ff_hard_%d'%nH] = np.array(fluxfactor_h)
        tbl['ff_all_%d'%nH] = np.array(fluxfactor_a)

        tbl['fr_soft_%d'%nH] = np.array(ratefactor_s)
        tbl['fr_hard_%d'%nH] = np.array(ratefactor_h)
        tbl['fr_all_%d'%nH] = np.array(ratefactor_a)

    tbl.write(outname, overwrite=True)


inname = '/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/narrowlineLRD_cat_mwnh.fits'
outname = '/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/narrowlineLRD_cat_mwnh_ff.fits'


get_convert_factor_sample(inname, outname)
