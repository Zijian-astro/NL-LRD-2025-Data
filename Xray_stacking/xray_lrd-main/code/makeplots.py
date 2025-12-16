import os, sys

datadir = '/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/'
sys.path.append(datadir+'output/')
from read_cstack_results import sourcerate

import numpy as np
from astropy.table import Table, join

import matplotlib.pyplot as plt

# plt.rcParams.update({
#     "text.usetex": True,
# #    "font.family": "Helvetica"
#     "font.family": "serif",
# })
plt.style.use('/Users/zijianzhang/Astro_Data/paper.mplstyle')


def read_and_match():
    origtbl = Table.read(datadir+'/narrowlineLRD_cat_good_nH_pimms.fits')
    # origtbl.remove_columns(['ff_soft', 'ff_hard', 'fr_soft', 'fr_hard', 'pimms_s', 'pimms_h'])

    finalgood = origtbl.copy()

    # read the fluxes
    tbl_soft_all_obs = Table.read(datadir+'/output/all_soft_results.fits')
    tbl_hard_all_obs = Table.read(datadir+'/output/all_hard_results.fits')

    exptimelist_soft = []
    exptimelist_hard = []
    softlimlist = []
    softv50list = []
    softv16list = []
    softv84list = []

    hardlimlist = []
    hardv50list = []
    hardv16list = []
    hardv84list = []

    alllimlist = []
    allv50list = []
    allv16list = []
    allv84list = []

    for index in range(len(finalgood)):
        name = finalgood['names'][index]
        print(name)
        subtbl_soft = tbl_soft_all_obs[tbl_soft_all_obs['names']==name]
        subtbl_hard = tbl_hard_all_obs[tbl_hard_all_obs['names']==name]

        exptimelist_soft.append(np.sum(subtbl_soft['expsrc']))
        exptimelist_hard.append(np.sum(subtbl_hard['expsrc']))

        rslist_s, prob_s, cdf_s = np.loadtxt(datadir+'/distributions/soft/%s.txt'%name)
        rslist_h, prob_h, cdf_h = np.loadtxt(datadir+'/distributions/hard/%s.txt'%name)
        # rslist_a, prob_a, cdf_a = np.loadtxt(datadir+'/distributions/all/%s.txt'%name)

        distribution_s = sourcerate(x=rslist_s,p=prob_s,c=cdf_s/np.max(cdf_s))
        distribution_h = sourcerate(x=rslist_h,p=prob_h,c=cdf_h/np.max(cdf_h))
        # distribution_a = sourcerate(x=rslist_a,p=prob_a,c=cdf_a/np.max(cdf_a))

        print(np.max(cdf_s), np.max(cdf_h)) #, np.max(cdf_a)

        softlimlist.append(distribution_s._ppf(0.9987))
        softv50list.append(distribution_s._ppf(0.50))
        softv84list.append(distribution_s._ppf(0.84))
        softv16list.append(distribution_s._ppf(0.16))
        hardlimlist.append(distribution_h._ppf(0.9987))
        hardv50list.append(distribution_h._ppf(0.50))
        hardv84list.append(distribution_h._ppf(0.84))
        hardv16list.append(distribution_h._ppf(0.16))
        # alllimlist.append(distribution_a._ppf(0.9987))
#         allv50list.append(distribution_a._ppf(0.50))
#         allv84list.append(distribution_a._ppf(0.84))
#         allv16list.append(distribution_a._ppf(0.16))

    finalgood['exptime_soft'] = exptimelist_soft
    finalgood['exptime_hard'] = exptimelist_hard

    finalgood['softlim'] = softlimlist
    finalgood['softv50'] = softv50list
    finalgood['softv84'] = softv84list
    finalgood['softv16'] = softv16list
    finalgood['hardlim'] = hardlimlist
    finalgood['hardv50'] = hardv50list
    finalgood['hardv84'] = hardv84list
    finalgood['hardv16'] = hardv16list
    # finalgood['alllim'] = alllimlist
#     finalgood['allv50'] = allv50list
#     finalgood['allv84'] = allv84list
#     finalgood['allv16'] = allv16list


    finalgood.write(datadir+'/all_lrds_final.fits', overwrite=True)

def print_stacked_rates():
    rslist_s, prob_s, cdf_s = np.loadtxt(datadir+'/distributions/soft/stacked2.txt')
    rslist_h, prob_h, cdf_h = np.loadtxt(datadir+'/distributions/hard/stacked2.txt')
    # rslist_a, prob_a, cdf_a = np.loadtxt(datadir+'/distributions/all/stacked2.txt')

    distribution_s = sourcerate(x=rslist_s,p=prob_s,c=cdf_s/np.max(cdf_s))
    distribution_h = sourcerate(x=rslist_h,p=prob_h,c=cdf_h/np.max(cdf_h))
    # distribution_a = sourcerate(x=rslist_a,p=prob_a,c=cdf_a/np.max(cdf_a))

    soft50, softel, softeu = \
            distribution_s._ppf(0.5), distribution_s._ppf(0.5)-distribution_s._ppf(0.16), distribution_s._ppf(0.84)-distribution_s._ppf(0.5)
    hard50, hardel, hardeu = \
            distribution_h._ppf(0.5), distribution_h._ppf(0.5)-distribution_h._ppf(0.16), distribution_h._ppf(0.84)-distribution_h._ppf(0.5)
    # all50, allel, alleu = \
 # distribution_a._ppf(0.5), distribution_a._ppf(0.5)-distribution_a._ppf(0.16), distribution_a._ppf(0.84)-distribution_a._ppf(0.5)

    print(r'%.2f_{-%.2f}^{+%.2f}'%(soft50*1e6, softel*1e6, softeu*1e6))
    print(r'%.2f_{-%.2f}^{+%.2f}'%(hard50*1e6, hardel*1e6, hardeu*1e6))
    # print(r'%.2f_{-%.2f}^{+%.2f}'%(all50*1e6, allel*1e6, alleu*1e6))


def get_expected_rates(nH):

    tbl_final = Table.read(datadir+'/all_lrds_final.fits')

    f_log = 0.83
    logLx = tbl_final['LHa']*f_log + 8.35
    logLHa_eu = tbl_final['LHa_eu']
    logLHa_el = tbl_final['LHa_el']

    Lx = 10**logLx

    softflux_exp = Lx * tbl_final['ff_soft_%d'%nH]
    hardflux_exp = Lx * tbl_final['ff_hard_%d'%nH]
    # allflux_exp =  Lx * tbl_final['ff_all_%d'%nH]

    # compute the expected error of expected fluxes
    softflux_exp_eu = softflux_exp * np.log(10) * logLHa_eu * f_log
    softflux_exp_el = softflux_exp * np.log(10) * logLHa_el * f_log

    hardflux_exp_eu = hardflux_exp * np.log(10) * logLHa_eu * f_log
    hardflux_exp_el = hardflux_exp * np.log(10) * logLHa_el * f_log

    # allflux_exp_eu = allflux_exp * np.log(10) * logLHa_eu * f_log
#     allflux_exp_el = allflux_exp * np.log(10) * logLHa_el * f_log

    softrate_exp = softflux_exp / tbl_final['pimms_s_%d'%nH]
    hardrate_exp = hardflux_exp / tbl_final['pimms_h_%d'%nH]
    # allrate_exp = allflux_exp / tbl_final['pimms_a_%d'%nH]

    softrate_exp_err = [softflux_exp_el / tbl_final['pimms_s_%d'%nH],\
                        softflux_exp_eu / tbl_final['pimms_s_%d'%nH]]
    hardrate_exp_err = [hardflux_exp_el / tbl_final['pimms_h_%d'%nH],\
                        hardflux_exp_eu / tbl_final['pimms_h_%d'%nH]]
    # allrate_exp_err = [allflux_exp_el / tbl_final['pimms_a_%d'%nH],\
    #                     allflux_exp_eu / tbl_final['pimms_a_%d'%nH]]

    uplims = np.ones(len(softrate_exp))

    weight = np.array(tbl_final['exptime']) / np.sum(tbl_final['exptime'])

    softrate_exp_mean = np.sum(softrate_exp * weight)
    hardrate_exp_mean = np.sum(hardrate_exp * weight)
    # allrate_exp_mean = np.sum(allrate_exp * weight)

    softrate_exp_eu = np.sqrt(np.sum(softflux_exp_eu**2*weight**2))
    softrate_exp_el = np.sqrt(np.sum(softflux_exp_el**2*weight**2))
    hardrate_exp_eu = np.sqrt(np.sum(hardflux_exp_eu**2*weight**2))
    hardrate_exp_el = np.sqrt(np.sum(hardflux_exp_el**2*weight**2))
    # allrate_exp_eu = np.sqrt(np.sum(allflux_exp_eu**2*weight**2))
    # allrate_exp_el = np.sqrt(np.sum(allflux_exp_el**2*weight**2))

    infodir = {'soft': [softrate_exp, softrate_exp_err, softrate_exp_mean, [[softrate_exp_eu], [softrate_exp_el]]],\
                'hard': [hardrate_exp, hardrate_exp_err, hardrate_exp_mean, [[hardrate_exp_eu], [hardrate_exp_el]]],\
                }
                #'all': [allrate_exp, allrate_exp_err, allrate_exp_mean, [[allrate_exp_eu], [allrate_exp_el]]]
    return infodir

# def plot_fluxes():
#     tbl_final = Table.read(datadir+'/all_lrds_final.fits')
#
#     infodir_nh21 = get_expected_rates(21)
#     infodir_nh22 = get_expected_rates(22)
#     infodir_nh23 = get_expected_rates(23)
#
#     # get the upper limit of the whole sample
#     rslist_s, prob_s, cdf_s = np.loadtxt(datadir+'/distributions/soft/stacked2.txt')
#     cdf_s = cdf_s / np.max(cdf_s)
#     dist_s = sourcerate(x=rslist_s, p=prob_s, c=cdf_s)
#     stack_uplim_s = dist_s._ppf(0.9987)
#
#     rslist_h, prob_h, cdf_h = np.loadtxt(datadir+'/distributions/hard/stacked2.txt')
#     cdf_h = cdf_h / np.max(cdf_h)
#     dist_h = sourcerate(x=rslist_h, p=prob_h, c=cdf_h)
#     stack_uplim_h = dist_h._ppf(0.9987)
#
#  #   rslist_a , prob_a, cdf_a = np.loadtxt(datadir+'/distributions/all/stacked2.txt')
# #     cdf_a = cdf_a / np.max(cdf_a)
# #     dist_a = sourcerate(x=rslist_a, p=prob_a, c=cdf_a)
# #     stack_uplim_a = dist_a._ppf(0.9987)
#
#     # plot a 3*3 plot
#     fig, axes = plt.subplots(3,3,figsize=[10,9])
#
#     # plot
#     yerr_plot_s = 0.2*tbl_final['softlim']
#     yerr_plot_h = 0.2*tbl_final['hardlim']
#     # yerr_plot_a = 0.2*tbl_final['alllim']
#
#     allinfolist = [infodir_nh21, infodir_nh22, infodir_nh23]
#
#     for infodir, axrow in zip(allinfolist, axes):
#         axrow[0].errorbar(infodir['soft'][0], tbl_final['softlim'],\
#                    xerr=infodir['soft'][1], yerr=yerr_plot_s,\
#                    uplims=[True]*len(tbl_final), fmt='o', color='k', label='LRDs', mfc='none')
#         axrow[0].errorbar(infodir['soft'][2],  stack_uplim_s, xerr=infodir['soft'][3],\
#                     yerr=[stack_uplim_s*0.2], color='r',fmt='s',\
#                    uplims=[True], label='Sample Average') # add upper limits
#         axes[0][0].set_title('0.5-2 keV', fontsize=16)
#
#         axrow[1].errorbar(infodir['hard'][0], tbl_final['hardlim'],\
#                    xerr=infodir['hard'][1], yerr=yerr_plot_h,\
#                    uplims=[True]*len(tbl_final), fmt='o', color='k', label='LRDs', mfc='none')
#         axrow[1].errorbar(infodir['hard'][2],  stack_uplim_h, xerr=infodir['hard'][3],\
#                     yerr=[stack_uplim_h*0.2], color='r',fmt='s',\
#                    uplims=[True], label='Sample Average') # add upper limits
#         axes[0][1].set_title('2-8 keV', fontsize=16)
#         #
#         # axrow[2].errorbar(infodir['all'][0], tbl_final['alllim'],\
#         #            xerr=infodir['all'][1], yerr=yerr_plot_a,\
#         #            uplims=[True]*len(tbl_final), fmt='o', color='k', label='LRDs', mfc='none')
#         # axrow[2].errorbar(infodir['all'][2],  stack_uplim_a, xerr=infodir['all'][3],\
#         #             yerr=[stack_uplim_a*0.2], color='r',fmt='s',\
#         #            uplims=[True], label='Sample Average') # add upper limits
#         # axes[0][2].set_title('0.5-8 keV', fontsize=16)
#
#     for axrow in axes:
#         for ax in axrow:
#             xlist = np.arange(3e-7, 3e-4, 3e-7)
#             ax.plot(xlist, xlist, 'k--')
#             ax.set_ylim([1e-7, 1e-4])
#             ax.set_xlim([10e-7, 10e-4])
#             ax.set_xscale('log')
#             ax.set_yscale('log')
#             ax.set_aspect('equal')
#
#         axrow[1].set_ylim([2e-7, 2e-4])
#         axrow[1].set_xlim([4e-7, 4e-4])
#
#     # labeling NH
#     nhlist = [21,22,23]
#     for nh, axrow in zip(nhlist, axes):
#         for ax in axrow:
#             ax.text(.01, .99, r'$\mathrm{N_H=10^{%d}cm^{-2}}$'%nh, ha='left', va='top', transform=ax.transAxes, fontsize=12)
#             ax.legend(frameon=True, fontsize=12, loc='lower left')
#             ax.tick_params(axis="both", direction="in", which="both", labelsize=12)
#
#     fig.supylabel(r'Observed Rate (count s$^{-1}$)', fontsize=16)
#     fig.supxlabel(r'Expected Rate (count s$^{-1}$)', fontsize=16)
#     plt.tight_layout()
#     plt.savefig('upperlim.pdf')
#     plt.show()


def plot_fluxes_only_nh22():
    tbl_final = Table.read(datadir + '/all_lrds_final.fits')

    infodir = get_expected_rates(22)  # Only NH=22

    # Stack upper limits
    rslist_s, prob_s, cdf_s = np.loadtxt(datadir + '/distributions/soft/stacked2.txt')
    cdf_s = cdf_s / np.max(cdf_s)
    dist_s = sourcerate(x=rslist_s, p=prob_s, c=cdf_s)
    stack_uplim_s = dist_s._ppf(0.9987)

    rslist_h, prob_h, cdf_h = np.loadtxt(datadir + '/distributions/hard/stacked2.txt')
    cdf_h = cdf_h / np.max(cdf_h)
    dist_h = sourcerate(x=rslist_h, p=prob_h, c=cdf_h)
    stack_uplim_h = dist_h._ppf(0.9987)

    fig, axes = plt.subplots(2, 1, figsize=[5, 9])

    yerr_plot_s = 0.2 * tbl_final['softlim']
    yerr_plot_h = 0.2 * tbl_final['hardlim']

    # Soft band (0.5–2 keV)
    ax_s = axes[0]
    ax_s.errorbar(infodir['soft'][0], tbl_final['softlim'],
                  xerr=infodir['soft'][1], yerr=yerr_plot_s,
                  uplims=[True]*len(tbl_final), fmt='o', color='k', label='Narrow-line LRDs', mfc='none')
    ax_s.errorbar(infodir['soft'][2], stack_uplim_s, xerr=infodir['soft'][3],
                  yerr=[stack_uplim_s * 0.2], color='r', fmt='s',
                  uplims=[True], label='Stacked')
    ax_s.set_title('0.5–2 keV', fontsize=15)

    # Hard band (2–8 keV)
    ax_h = axes[1]
    ax_h.errorbar(infodir['hard'][0], tbl_final['hardlim'],
                  xerr=infodir['hard'][1], yerr=yerr_plot_h,
                  uplims=[True]*len(tbl_final), fmt='o', color='k', label='Narrow-line LRDs', mfc='none')
    ax_h.errorbar(infodir['hard'][2], stack_uplim_h, xerr=infodir['hard'][3],
                  yerr=[stack_uplim_h * 0.2], color='r', fmt='s',
                  uplims=[True], label='Stacked')
    ax_h.set_title('2–8 keV', fontsize=15)

    for ax in axes:
        xlist = np.arange(3e-7, 3e-4, 3e-7)
        ax.plot(xlist, xlist, 'k--')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_aspect('equal')
        ax.tick_params(axis="both", direction="in", which="both", labelsize=17)
        ax.legend(frameon=True, fontsize=15, loc='lower right')
        ax.text(0.01, 0.99, r'$N_\mathrm{H}=10^{22}\,\mathrm{cm}^{-2}$', 
                transform=ax.transAxes, fontsize=12, va='top', ha='left')

    axes[0].set_ylim([1e-7, 1e-4])
    axes[0].set_xlim([1e-6, 1e-3])
    axes[1].set_ylim([2e-7, 2e-4])
    axes[1].set_xlim([4e-7, 4e-4])
    
    axes[0].text(0.88, 0.93, '(a)', fontsize=15, transform=axes[0].transAxes, color='black', weight='bold')
    axes[1].text(0.88, 0.93, '(b)', fontsize=15, transform=axes[1].transAxes, color='black', weight='bold')
    axes[1].set_xlabel(r'Expected Rate (count s$^{-1}$)', fontsize=15)
    fig.supylabel(r'Observed Rate (count s$^{-1}$)', fontsize=15, x=0.03)
    fig.subplots_adjust(hspace=0.2)
    # plt.tight_layout()
    plt.savefig('upperlim_nh22_vertical.pdf',  bbox_inches='tight')
    plt.show()

def read_jin12():
    tbl1 = Table.read(datadir+'/../jin12_table1.txt', format='ascii')
    logHa = np.loadtxt(datadir+'/../jin12.txt')
    logHa_broad = logHa[:,0]
    logHa_narrow = logHa[:,1]


    tbl1['logLx'] = np.log10(tbl1['col4'])+44
    tbl1['logHa'] = logHa_broad#np.log10(10**logHa_N+10**logHa_Y)

    return tbl1

def plot_LxLHa():
    # tbl_jin12 = read_jin12()

    # read the LRDs
    tbl_final = Table.read(datadir+'/all_lrds_final.fits')
    tbl = tbl_final

    logLx = 1.11*tbl_final['LHa']-3.5
    logLHa_eu = tbl_final['LHa_eu']
    logLHa_el = tbl_final['LHa_el']

    Lx = 10**logLx

    # convert count limits to luminosity limits
    soft_lum_lim_21 = tbl['softlim'] * tbl['pimms_s_21'] / tbl['ff_soft_21']
    hard_lum_lim_21 = tbl['hardlim'] * tbl['pimms_h_21'] / tbl['ff_hard_21']
    all_lum_lim_21 = tbl['hardlim'] * tbl['pimms_h_21'] / tbl['ff_hard_21']
    lumlim_21 = np.min([soft_lum_lim_21, hard_lum_lim_21, all_lum_lim_21], axis=0)

    soft_lum_lim_22 = tbl['softlim'] * tbl['pimms_s_22'] / tbl['ff_soft_22']
    hard_lum_lim_22 = tbl['hardlim'] * tbl['pimms_h_22'] / tbl['ff_hard_22']
    all_lum_lim_22 = tbl['hardlim'] * tbl['pimms_h_22'] / tbl['ff_hard_22']
    lumlim_22 = np.min([soft_lum_lim_22, hard_lum_lim_22, all_lum_lim_22], axis=0)

    soft_lum_lim_23 = tbl['softlim'] * tbl['pimms_s_23'] / tbl['ff_soft_23']
    hard_lum_lim_23 = tbl['hardlim'] * tbl['pimms_h_23'] / tbl['ff_hard_23']
    all_lum_lim_23 = tbl['hardlim'] * tbl['pimms_h_23'] / tbl['ff_hard_23']
    lumlim_23 = np.min([soft_lum_lim_23, hard_lum_lim_23, all_lum_lim_23], axis=0)

    logLxlim_21 = np.log10(lumlim_21)
    logLxlim_22 = np.log10(lumlim_22)
    logLxlim_23 = np.log10(lumlim_23)


    LHa_arr = np.arange(39, 45, 0.1)

    # Jin+12 relation
    Lx_arr = 0.83 * LHa_arr + 8.35
    yerr_plot = [0.2] * len(logLxlim_21)
    uplims = [True] * len(logLxlim_21)

    fig, axes = plt.subplots(ncols=3,figsize=[10,3.7])

    for ax in axes:
        # ax.plot(tbl_jin12['logHa'], tbl_jin12['logLx'], '+', color='gray',\
      #          label='Low-z AGNs (Jin+12)')
        ax.plot(LHa_arr, Lx_arr, 'k--')

        ax.set_xlim([41, 44.5])
        ax.set_ylim([42, 45.5])

        ax.tick_params(axis="both", direction="in", which="both", labelsize=12)

    axes[0].errorbar(tbl_final['LHa'], logLxlim_21,\
               xerr=[logLHa_eu, logLHa_el], yerr=yerr_plot,\
               uplims=uplims, fmt='o', color='k', label='LRDs (Yue+24)',\
                mfc='none')
    axes[1].errorbar(tbl_final['LHa'], logLxlim_22,\
               xerr=[logLHa_eu, logLHa_el], yerr=yerr_plot,\
               uplims=uplims, fmt='o', color='k', label='LRDs (Yue+24)',\
                mfc='none')
    axes[2].errorbar(tbl_final['LHa'], logLxlim_23,\
               xerr=[logLHa_eu, logLHa_el], yerr=yerr_plot,\
               uplims=uplims, fmt='o', color='k', label='LRDs (Yue+24)',\
                mfc='none')
    axes[1].set_xlabel(r'$L_{\mathrm{H}{\alpha}}$ (erg s$^{-1}$)', fontsize=14)
    axes[0].set_ylabel(r'$L_{X}[2-10\mathrm{keV}]$ (erg s$^{-1}$)', fontsize=14)

    axes[0].legend(fontsize=12, frameon=False)

    axes[0].set_title(r'$\mathrm{N_H=10^{21}cm^{-2}}$', fontsize=12)
    axes[1].set_title(r'$\mathrm{N_H=10^{22}cm^{-2}}$', fontsize=12)
    axes[2].set_title(r'$\mathrm{N_H=10^{23}cm^{-2}}$', fontsize=12)

    plt.tight_layout()
    plt.savefig('LxLHa.pdf')
    plt.show()

def main():
    read_and_match()
    plot_fluxes_only_nh22()
    # plot_LxLHa()

if __name__=='__main__':
    main()
