import os, sys
import numpy as np
from astropy.table import Table

def generate_latex_table_code(filename):
    tbl = Table.read(filename)

    # tbl1: columns are name, RA, Dec, redshift, logLHa, FWHM, logMbh, ref

    for index in range(len(tbl)):
        name = tbl['names'][index]
        name = name.replace('_', '-')

        ra = tbl['RA'][index]
        dec = tbl['Dec'][index]
        redshift = tbl['redshift'][index]

        logLHa = tbl['LHa'][index]
        if logLHa>50:
            logLHa = logLHa-42
        logLHa_ue = tbl['LHa_eu'][index]
        logLHa_le = tbl['LHa_el'][index]

        FWHM = tbl['FWHMHa'][index]
        FWHMeu = tbl['FWHMHa_eu'][index]
        FWHMel = tbl['FWHMHa_el'][index]

        BHmass = tbl['BHmass'][index]
        BHmass_eu = tbl['BHmass_eu'][index]
        BHmass_el = tbl['BHmass_el'][index]


        string = r'%s & %.5f & %.5f & %.3f & $%.2f^{+%.2f}_{-%.2f}$ & $%d^{+%d}_{-%d}$ & $%.2f^{+%.2f}_{-%.2f}$\\'%(name, ra, dec, redshift, logLHa, logLHa_ue, logLHa_le, FWHM, FWHMeu, FWHMel, BHmass, BHmass_eu, BHmass_el)
        print(string)

def generate_latex_table_2_code(filename):
    tbl = Table.read(filename)

    # tbl2: columns are name, cts_soft, cts_hard, flux_soft, flux_hard, lum_soft, lum_hard

    print(tbl.colnames)

    for index in range(len(tbl)):
        name = tbl['names'][index]
        name = name.replace('_', '-')

        soft_count = tbl['softv50'][index] / 0.9
        soft_count_eu = (tbl['softv84'][index] - soft_count)/0.9
        soft_count_el = (soft_count - tbl['softv16'][index])/0.9

        hard_count = tbl['hardv50'][index] / 0.9
        hard_count_eu = (tbl['hardv84'][index] - hard_count) / 0.9
        hard_count_el = (-tbl['hardv16'][index] + hard_count) / 0.9

        soft_flux = soft_count * tbl['pimms_s'][index]
        soft_flux_eu = soft_count_eu * tbl['pimms_s'][index]
        soft_flux_el = soft_count_el * tbl['pimms_s'][index]

        hard_flux = hard_count * tbl['pimms_h'][index]
        hard_flux_eu = hard_count_eu * tbl['pimms_h'][index]
        hard_flux_el = hard_count_el * tbl['pimms_h'][index]

        soft_lum = soft_flux / tbl['ff_soft'][index]
        soft_lum_eu = soft_flux_eu / tbl['ff_soft'][index]
        soft_lum_el = soft_flux_el / tbl['ff_soft'][index]

        hard_lum = hard_flux / tbl['ff_hard'][index]
        hard_lum_eu = hard_flux_eu / tbl['ff_hard'][index]
        hard_lum_el = hard_flux_el / tbl['ff_hard'][index]

        soft_lum_lim = tbl['softlim'][index] * tbl['pimms_s'][index] / tbl['ff_soft'][index]
        hard_lum_lim = tbl['hardlim'][index] * tbl['pimms_h'][index] / tbl['ff_hard'][index]

        line = r'%s & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $<%.2f$ & $<%.2f$ \\'\
                %(name, soft_count*1e6, soft_count_eu*1e6, soft_count_el*1e6,\
                  hard_count*1e6, hard_count_eu*1e6, hard_count_el*1e6,\
                  soft_flux*1e17, soft_flux_eu*1e17, soft_flux_el*1e17,\
                  hard_flux*1e17, hard_flux_eu*1e17, hard_flux_el*1e17,\
                soft_lum_lim*1e-43, hard_lum_lim*1e-43)
#                  soft_lum/1e42, soft_lum_eu/1e42, soft_lum_el*1e42,\
#                  hard_lum/1e42, hard_lum_eu/1e42,, hard_lum_el*1e42)

        print(line)

def generate_table_2_new(filename):
    tbl = Table.read(filename)

    # tbl2: columns are name, cts_soft, cts_hard, flux_soft, flux_hard, flux_all, lum_soft, lum_hard, lum_all
    print(tbl.colnames)

    for index in range(len(tbl)):
        name = tbl['names'][index]
        name = name.replace('_', '-')

        ra = tbl['RA'][index]
        dec = tbl['Dec'][index]

        if ra>200:#AEGIS
            field='AEGIS'
        elif ra>150:#CDFN
            field='CDF-N'
        else:# CDFS
            field='CDF-S'

        exptime = tbl['exptime_soft'][index]
        nhmw = tbl['NH_MW'][index]

        soft_count = tbl['softv50'][index]
        soft_count_eu = (tbl['softv84'][index] - soft_count)
        soft_count_el = (soft_count - tbl['softv16'][index])

        hard_count = tbl['hardv50'][index]
        hard_count_eu = (tbl['hardv84'][index] - hard_count)
        hard_count_el = (-tbl['hardv16'][index] + hard_count)

        all_count = tbl['allv50'][index]
        all_count_eu = (tbl['allv84'][index] - all_count)
        all_count_el = (-tbl['allv16'][index] + all_count)

        line = r'{} & {} & {:.2f} & {:.2f} &'.format(name, field, nhmw/1e20, exptime/1e6)
        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & '.format(\
                    soft_count*1e6, soft_count_eu*1e6, soft_count_el*1e6)
        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ &'.format(\
                hard_count*1e6, hard_count_eu*1e6, hard_count_el*1e6)
        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ \\'.format(\
                all_count*1e6, all_count_eu*1e6, all_count_el*1e6)

        print(line)

    print(np.sum(tbl['exptime_soft']))


def generate_table_3_new(filename):
    tbl = Table.read(filename)

    # tbl2: columns are name, cts_soft, cts_hard, flux_soft, flux_hard, flux_all, lum_soft, lum_hard, lum_all
    print(tbl.colnames)


    for index in range(len(tbl)):
        name = tbl['names'][index]
        name = name.replace('_', '-')

        soft_count = tbl['softv50'][index]
        soft_count_eu = (tbl['softv84'][index] - soft_count)
        soft_count_el = (soft_count - tbl['softv16'][index])

        hard_count = tbl['hardv50'][index]
        hard_count_eu = (tbl['hardv84'][index] - hard_count)
        hard_count_el = (-tbl['hardv16'][index] + hard_count)

        all_count = tbl['allv50'][index]
        all_count_eu = (tbl['allv84'][index] - all_count)
        all_count_el = (-tbl['allv16'][index] + all_count)

        soft_flux_21 = soft_count * tbl['pimms_s_21'][index]
        soft_flux_eu_21 = soft_count_eu * tbl['pimms_s_21'][index]
        soft_flux_el_21 = soft_count_el * tbl['pimms_s_21'][index]

        all_flux_21 = all_count * tbl['pimms_h_21'][index]
        all_flux_eu_21 = all_count_eu * tbl['pimms_h_21'][index]
        all_flux_el_21 = all_count_el * tbl['pimms_h_21'][index]

        hard_flux_21 = hard_count * tbl['pimms_h_21'][index]
        hard_flux_eu_21 = hard_count_eu * tbl['pimms_h_21'][index]
        hard_flux_el_21 = hard_count_el * tbl['pimms_h_21'][index]

        soft_lum_21 = soft_flux_21 / tbl['ff_soft_21'][index]
        soft_lum_eu_21 = soft_flux_eu_21 / tbl['ff_soft_21'][index]
        soft_lum_el_21 = soft_flux_el_21 / tbl['ff_soft_21'][index]

        hard_lum_21 = hard_flux_21 / tbl['ff_hard_21'][index]
        hard_lum_eu_21 = hard_flux_eu_21 / tbl['ff_hard_21'][index]
        hard_lum_el_21 = hard_flux_el_21 / tbl['ff_hard_21'][index]

        all_lum_21 = all_flux_21 / tbl['ff_all_21'][index]
        all_lum_eu_21 = all_flux_eu_21 / tbl['ff_all_21'][index]
        all_lum_el_21 = all_flux_el_21 / tbl['ff_all_21'][index]

        soft_lum_lim_21 = tbl['softlim'][index] * tbl['pimms_s_21'][index] / tbl['ff_soft_21'][index]
        hard_lum_lim_21 = tbl['hardlim'][index] * tbl['pimms_h_21'][index] / tbl['ff_hard_21'][index]
        all_lum_lim_21 = tbl['alllim'][index] * tbl['pimms_h_21'][index] / tbl['ff_hard_21'][index]

        soft_flux_22 = soft_count * tbl['pimms_s_22'][index]
        soft_flux_eu_22 = soft_count_eu * tbl['pimms_s_22'][index]
        soft_flux_el_22 = soft_count_el * tbl['pimms_s_22'][index]

        all_flux_22 = all_count * tbl['pimms_h_22'][index]
        all_flux_eu_22 = all_count_eu * tbl['pimms_h_22'][index]
        all_flux_el_22 = all_count_el * tbl['pimms_h_22'][index]

        hard_flux_22 = hard_count * tbl['pimms_h_22'][index]
        hard_flux_eu_22 = hard_count_eu * tbl['pimms_h_22'][index]
        hard_flux_el_22 = hard_count_el * tbl['pimms_h_22'][index]

        soft_lum_22 = soft_flux_22 / tbl['ff_soft_22'][index]
        soft_lum_eu_22 = soft_flux_eu_22 / tbl['ff_soft_22'][index]
        soft_lum_el_22 = soft_flux_el_22 / tbl['ff_soft_22'][index]

        hard_lum_22 = hard_flux_22 / tbl['ff_hard_22'][index]
        hard_lum_eu_22 = hard_flux_eu_22 / tbl['ff_hard_22'][index]
        hard_lum_el_22 = hard_flux_el_22 / tbl['ff_hard_22'][index]

        all_lum_22 = all_flux_22 / tbl['ff_all_22'][index]
        all_lum_eu_22 = all_flux_eu_22 / tbl['ff_all_22'][index]
        all_lum_el_22 = all_flux_el_22 / tbl['ff_all_22'][index]

        soft_lum_lim_22 = tbl['softlim'][index] * tbl['pimms_s_22'][index] / tbl['ff_soft_22'][index]
        hard_lum_lim_22 = tbl['hardlim'][index] * tbl['pimms_h_22'][index] / tbl['ff_hard_22'][index]
        all_lum_lim_22 = tbl['alllim'][index] * tbl['pimms_h_22'][index] / tbl['ff_hard_22'][index]

        soft_flux_23 = soft_count * tbl['pimms_s_23'][index]
        soft_flux_eu_23 = soft_count_eu * tbl['pimms_s_23'][index]
        soft_flux_el_23 = soft_count_el * tbl['pimms_s_23'][index]

        all_flux_23 = all_count * tbl['pimms_h_23'][index]
        all_flux_eu_23 = all_count_eu * tbl['pimms_h_23'][index]
        all_flux_el_23 = all_count_el * tbl['pimms_h_23'][index]

        hard_flux_23 = hard_count * tbl['pimms_h_23'][index]
        hard_flux_eu_23 = hard_count_eu * tbl['pimms_h_23'][index]
        hard_flux_el_23 = hard_count_el * tbl['pimms_h_23'][index]

        soft_lum_23 = soft_flux_23 / tbl['ff_soft_23'][index]
        soft_lum_eu_23 = soft_flux_eu_23 / tbl['ff_soft_23'][index]
        soft_lum_el_23 = soft_flux_el_23 / tbl['ff_soft_23'][index]

        hard_lum_23 = hard_flux_23 / tbl['ff_hard_23'][index]
        hard_lum_eu_23 = hard_flux_eu_23 / tbl['ff_hard_23'][index]
        hard_lum_el_23 = hard_flux_el_23 / tbl['ff_hard_23'][index]

        all_lum_23 = all_flux_23 / tbl['ff_all_23'][index]
        all_lum_eu_23 = all_flux_eu_23 / tbl['ff_all_23'][index]
        all_lum_el_23 = all_flux_el_23 / tbl['ff_all_23'][index]

        soft_lum_lim_21 = tbl['softlim'][index] * tbl['pimms_s_21'][index] / tbl['ff_soft_21'][index]
        hard_lum_lim_21 = tbl['hardlim'][index] * tbl['pimms_h_21'][index] / tbl['ff_hard_21'][index]
        all_lum_lim_21 = tbl['alllim'][index] * tbl['pimms_h_21'][index] / tbl['ff_hard_21'][index]
        lumlim_21 = np.min([soft_lum_lim_21, hard_lum_lim_21, all_lum_lim_21], axis=0)

        soft_lum_lim_22 = tbl['softlim'][index] * tbl['pimms_s_22'][index] / tbl['ff_soft_22'][index]
        hard_lum_lim_22 = tbl['hardlim'][index] * tbl['pimms_h_22'][index] / tbl['ff_hard_22'][index]
        all_lum_lim_22 = tbl['alllim'][index] * tbl['pimms_h_22'][index] / tbl['ff_hard_22'][index]
        lumlim_22 = np.min([soft_lum_lim_22, hard_lum_lim_22, all_lum_lim_22], axis=0)

        soft_lum_lim_23 = tbl['softlim'][index] * tbl['pimms_s_23'][index] / tbl['ff_soft_23'][index]
        hard_lum_lim_23 = tbl['hardlim'][index] * tbl['pimms_h_23'][index] / tbl['ff_hard_23'][index]
        all_lum_lim_23 = tbl['alllim'][index] * tbl['pimms_h_23'][index] / tbl['ff_hard_23'][index]
        lumlim_23 = np.min([soft_lum_lim_23, hard_lum_lim_23, all_lum_lim_23], axis=0)

        line = r'{} & '.format(name)
        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & $<{:.2f}$ & '.format(\
                soft_flux_21*1e17, soft_flux_eu_21*1e17, soft_flux_el_21*1e17,\
                hard_flux_21*1e17, hard_flux_eu_21*1e17, hard_flux_el_21*1e17,\
                all_flux_21*1e17, all_flux_eu_21*1e17, all_flux_el_21*1e17, lumlim_21*1e-43)

        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & $<{:.2f}$ & '.format(\
                soft_flux_22*1e17, soft_flux_eu_22*1e17, soft_flux_el_22*1e17,\
                hard_flux_22*1e17, hard_flux_eu_22*1e17, hard_flux_el_22*1e17,\
                all_flux_22*1e17, all_flux_eu_22*1e17, all_flux_el_22*1e17, lumlim_22*1e-43)

        line+= r'${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ & $<{:.2f}$\\'.format(\
                soft_flux_23*1e17, soft_flux_eu_23*1e17, soft_flux_el_23*1e17,\
                hard_flux_23*1e17, hard_flux_eu_23*1e17, hard_flux_el_23*1e17,\
                all_flux_23*1e17, all_flux_eu_23*1e17, all_flux_el_23*1e17, lumlim_23*1e-43)


        print(line)
        #& $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $%.2f^{+%.2f}_{-%.2f}$ & $<%.2f$ & $<%.2f$ \\'\
        #        %(name, soft_count*1e6, soft_count_eu*1e6, soft_count_el*1e6,\
        #          hard_count*1e6, hard_count_eu*1e6, hard_count_el*1e6,\
        #          soft_flux*1e17, soft_flux_eu*1e17, soft_flux_el*1e17,\
        #          hard_flux*1e17, hard_flux_eu*1e17, hard_flux_el*1e17,\
        #        soft_lum_lim*1e-43, hard_lum_lim*1e-43)



if __name__=='__main__':
#    generate_latex_table_code('./all_lrds_cstack_good.fits')
    generate_table_2_new('/Users/minghao/Research/Projects/JWST/LRDs/data/new/all_lrds_final.fits')


