import os, sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table, vstack
from scipy.stats import poisson
from scipy import stats

from scipy.special import gamma, gammaincc
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from scipy.signal import convolve
import glob

# edit this path when you try to rerun the code
dir_cstack_output = '/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/output/'
dir_distribution = os.path.abspath(dir_cstack_output+'/../distributions/')

'''
Utility functions
'''

def read_results(filename):
    '''
    Extracting information from CSTACK output file
    '''
    f = open(filename, 'r')

    nline = 0

    cntsrclist = []
    expsrclist = []
    pixsrclist = []

    cntbkglist = []
    expbkglist = []
    pixbkglist = []

    namelist = []
    bandlist = []
    nobslist = []

    radiuslist = []
    offanglelist = []
    bandlist = []

    while True:
        line = f.readline()

        valuelist = line.split()
        if nline==1:
            band = valuelist[6]
        else:
            try:
                nobslist.append(int(valuelist[1]))
                namelist.append(valuelist[0])
                cntsrclist.append(float(valuelist[3]))
                expsrclist.append(float(valuelist[4]))
                pixsrclist.append(float(valuelist[5]))
                cntbkglist.append(float(valuelist[6]))
                expbkglist.append(float(valuelist[7]))
                pixbkglist.append(float(valuelist[8]))

                radiuslist.append(float(valuelist[2]))
                offanglelist.append(float(valuelist[12]))
                bandlist.append(band)
                tbl = Table({'names': namelist,\
                             'nobs': nobslist, \
                             'cntsrc': cntsrclist,\
                             'expsrc': expsrclist,\
                             'pixsrc': pixsrclist,\
                             'cntbkg': cntbkglist,\
                             'expbkg': expbkglist,\
                             'pixbkg': pixbkglist,\
                             'band': bandlist})
            except:
                pass

        if valuelist[0]=='sum_all':
            break
        if not line:
            break

        nline +=1

    tbl = Table({'names': namelist,\
                 'nobs': nobslist,\
                 'rad': radiuslist,\
                 'offax': offanglelist,\
                 'cntsrc': cntsrclist,\
                 'expsrc': expsrclist,\
                 'pixsrc': pixsrclist,\
                 'cntbkg': cntbkglist,\
                 'expbkg': expbkglist,\
                 'pixbkg': pixbkglist,\
                 'band': bandlist})

    return tbl


def find_good_dynamical_range(rslist, problist, thres=0.01):
    '''
    return a proper range of rs given its probability distribution
    '''
    mask = problist>thres
    rslist_tr = rslist[mask]
    return np.max(np.abs([rslist_tr[0], rslist_tr[-1]]))

'''
The following functions is the core of computing count rate distribution.
'''
def prob_rs_list_obs(rslist, cntsrclist, expsrclist, pixsrclist,\
                     cntbkglist, expbkglist, pixbkglist, ecf):
    '''
    Compute the distribution of rs (for a list of rs)
    given a list of observations.
    '''
    # weighted by expsrclist

    # on this grid
    xgrid = rslist

    # this correspond to weight = expsrc
    Ttot = np.sum(expsrclist)
    Nstot = np.sum(cntsrclist)

    # prob of Ns/T
    # this part compute the probability distribution of r_all in Eq. 2
    coeff0 = 1/Ttot/ecf # coeff0 = rate/count, count = rate/coeff
    P0 = stats.gamma.pdf(xgrid/coeff0, Nstot+1) /coeff0
    Pall = P0

    # the following lines takes into account the background counts.
    for index in range(len(cntsrclist)):
        coeff = expsrclist[index] * pixsrclist[index] / expbkglist[index] / pixbkglist[index] / Ttot * (-1)/ecf
        Pi = stats.gamma.pdf(xgrid/coeff, cntbkglist[index]+1) / coeff * (-1)
        Pall = convolve(Pall, Pi, mode='same') * (rslist[1]-rslist[0])

    return Pall

'''
We define a new class to hold the derived distribution
'''
class sourcerate(stats.rv_continuous):
    def __init__(self, x, p, c=None):
        '''
        x is the grid of the variable (in this work, the count rate)
        p is the probabilistic distribution
        c is the cumulative distribution
        '''
        super().__init__(a=x.min(), b=x.max())
        self.x = x
        self.p = p

        if c is None:
            dx = x[1] - x[0]
            cdflist = [np.sum(self.p[:index])*dx for index in range(len(x))]
            self.c = np.array(cdflist)
        else:
            self.c = c

    def _cdf(self, x):
        return interp1d(self.x, self.c)(x)
    def _pdf(self, x):
        return interp1d(self.x,self.p)(x)
    def _ppf(self, q):
        return interp1d(self.c, self.x)(q)

def summarize_info(direct=dir_cstack_output):
    '''
    Summarize the output from cstack into tables for further analysis
    '''
    files_soft = glob.glob(direct+'/*/stat_cts_500_2000.out')
    files_soft.sort()
    tbllist_soft = [read_results(f) for f in files_soft]
    tbl_all_soft = vstack(tbllist_soft)
    tbl_all_soft.write(direct+'/all_soft_results.fits', overwrite=True)

    files_hard = glob.glob(direct+'/*/stat_cts_2000_8000.out')
    files_hard.sort()
    tbllist_hard = [read_results(f) for f in files_hard]
    tbl_all_hard = vstack(tbllist_hard)
    tbl_all_hard.write(direct+'/all_hard_results.fits', overwrite=True)


def save_all_prob_dists(tbl, savedir):
    '''
    Save the count rate distribution of all objects.
    Input:
        tblname: the saved information directory from summarize_info()
        savedir: where to save the probabilistic functions
    Output: None
    '''

    allnames = set(tbl['names'])

    for name in allnames:
        subtbl = tbl[tbl['names']==name]
        cntsrc, expsrc, pixsrc, cntbkg, expbkg, pixbkg = \
            subtbl['cntsrc'], subtbl['expsrc'], subtbl['pixsrc'],\
            subtbl['cntbkg'], subtbl['expbkg'], subtbl['pixbkg']

        '''
        Given the different background counts and source counts, the dynamical range
        of the derived source rate can be very different. The following lines ensures
        a good coverage of rs while keep the code fast enough.
        '''

        factor_soft = 1
        if name in ['CEERS_00672', 'CEERS_00717']:
            factor_soft *= 0.2

        elif name=='M23_10013704_2':
            factor_soft *= 0.2

        rslist = np.arange(-5e-4, 5e-4+1e-9, 1e-9)*factor_soft

        '''
        We now compute and save the probability distribution
        '''
        problist = prob_rs_list_obs(rslist, cntsrc, expsrc, pixsrc, cntbkg, expbkg, pixbkg, ecf=0.9)

        # also save the cdf for future use
        # to speed up, use a coarse grid

        xedge = find_good_dynamical_range(rslist, problist)
        xlist_coarse =  np.linspace(-xedge, xedge, 10000)
        pdf_coarse = np.interp(xlist_coarse, rslist, problist)

        cdf_coarse = [np.trapz(pdf_coarse[:index], xlist_coarse[:index])\
               for index in range(len(xlist_coarse))]

        np.savetxt(savedir+'/'+name+'.txt', [xlist_coarse, np.array(pdf_coarse), np.array(cdf_coarse)])

def print_all_prob_dist_limits(band='soft'):
    tbl_cstack = Table.read(dir_cstack_output+'/../all_lrds_final.fits')

    for index in range(len(tbl_cstack)):
        name = tbl_cstack['names'][index]

        dist_file = dir_cstack_output+'/../distributions/%s/%s.txt'%(band,name)
        rslist, problist, cdflist = np.loadtxt(dist_file)

        rsdist = sourcerate(x=rslist, p=problist, c=cdflist)

        try:
            low = rsdist._ppf(0.0013)
            med = rsdist._ppf(0.50)
            high = rsdist._ppf(0.9987)
            print(name, low, med, high)
        except:
            print(name, 'bad')
            plt.plot(rslist, rsdist.c)
            plt.show()

def distribution_of_sum(xlist, plist, weight=None):
    # find the correct dynamical range
    xedges = []
    dxlist = []

    if weight is None:
        weight= np.ones(len(xlist))
    weight = weight/np.sum(weight)

    xlist = [xlist[index]*weight[index] for index in range(len(xlist))]
    plist = [plist[index]/weight[index] for index in range(len(xlist))]

    for i, (x,p) in enumerate(zip(xlist,plist)):
        xedges.append(find_good_dynamical_range(x, p))
        dxlist.append(x[1]-x[0])

    xedge_master = np.max(xedges)*5
    Ngrid = int(2*xedge_master/np.min(dxlist)) + 1

    xlist_master = np.linspace(-xedge_master, xedge_master, Ngrid)
    dx_master = xlist_master[1]-xlist_master[0]

    plist_resample = [np.interp(xlist_master, x, p) for x, p in zip(xlist, plist)]
    psum = plist_resample[0]
    for index, p in enumerate(plist_resample[1:]):
        psum = convolve(p, psum, mode='same') * dx_master

    # this is too slow. We may want to do a coarser grid

    xlist_coarse =  np.linspace(-xedge_master, xedge_master, 10000)
    dx_coarse = xlist_master[1]-xlist_master[0]
    psum_coarse = np.interp(xlist_coarse, xlist_master, psum)

    cdf_coarse = [np.trapz(psum_coarse[:index], xlist_coarse[:index])\
               for index in range(len(xlist_coarse))]

    return xlist_coarse, psum_coarse, cdf_coarse

def get_master_distribution_new(band):
    # use exptime as the weight
    tbl_all = Table.read(dir_cstack_output+'/../narrowlineLRD_cat.fits')

    rslist_common = np.arange(-1e-4, 1e-4, 5e-9)
    problist_all = []
    rlist_all = []

    for index in range(len(tbl_all)):
        name = tbl_all['uid'][index]
#        read the distribution 
        rslist, problist, cdflist = np.loadtxt(dir_distribution+'/%s/%s.txt'%(band,name))
        rlist_all.append(rslist)
        problist_all.append(problist)

    rs_master, p_master, cdf_master = distribution_of_sum(rlist_all, problist_all, weight=tbl_all['exptime'])
    np.savetxt(dir_distribution+'/%s/stacked2.txt'%band,\
               [rs_master, np.array(p_master), np.array(cdf_master)])

    plt.plot(rs_master, p_master)
    plt.show()

def full_band_distribution(soft_dist_dir, hard_dist_dir):
    '''
    Get the distribution of full (soft+hard) photon rates.
    Simply computed as rs_full = rs_soft+rs_hard
    Hence convolving P_soft(rs) and P_hard(rs)
    '''

    allnames = os.listdir(soft_dist_dir)
    allnames.sort()

    for name in allnames:
        if not name.endswith('.txt'):
            continue
        elif name.startswith('full'):
            continue
        else:
            objname = name[:-4]
            rslist1, problist1, cdf1 = np.loadtxt(soft_dist_dir+name)
            rslist2, problist2, cdf2 = np.loadtxt(hard_dist_dir+name)

            redge1 = find_good_dynamical_range(rslist1, problist1)
            redge2 = find_good_dynamical_range(rslist2, problist2)

#            need to sample it to a common grid
            drs = np.min([rslist1[1]-rslist1[0], rslist2[1]-rslist2[0]])
            redge = np.max([redge1, redge2])*2

            rslist = np.linspace(-redge, redge, int(2*redge/drs))
            problist1_resample = np.interp(rslist, rslist1, problist1)
            problist2_resample = np.interp(rslist, rslist2, problist2)

            problist = convolve(problist1_resample, problist2_resample, mode='same')\
                    * (rslist[1]-rslist[0])

            cdflist = [np.trapz(problist[:index], rslist[:index]) for index in range(len(rslist))]

            np.savetxt(dir_cstack_output+'/../distributions/full/%s'%name,\
                       [rslist, np.array(problist), np.array(cdflist)])


def full_band_distribution_master():
    # simply soft + hard

    soft_dist_dir = dir_cstack_output+'/../distributions/soft/'
    hard_dist_dir = dir_cstack_output+'/../distributions/hard/'

    name = 'stacked.txt'
    rslist1, problist1, cdf1 = np.loadtxt(soft_dist_dir+name)
    rslist2, problist2, cdf2 = np.loadtxt(hard_dist_dir+name)

    redge1 = find_good_dynamical_range(rslist1, problist1)
    redge2 = find_good_dynamical_range(rslist2, problist2)

    drs = np.min([rslist1[1]-rslist1[0], rslist2[1]-rslist2[0]])
    redge = np.max([redge1, redge2])*2

    rslist = np.linspace(-redge, redge, int(2*redge/drs))
    problist1_resample = np.interp(rslist, rslist1, problist1)
    problist2_resample = np.interp(rslist, rslist2, problist2)

    problist = convolve(problist1_resample, problist2_resample, mode='same')\
            * (rslist[1]-rslist[0])

    cdflist = [np.trapz(problist[:index], rslist[:index]) for index in range(len(rslist))]

    np.savetxt(dir_cstack_output+'/../distributions/full/%s'%name,\
               [rslist, np.array(problist), np.array(cdflist)])


if __name__=='__main__':

    # first save individual CSTACK outputs to a master table
    summarize_info()

    # We can now read the soft and hard band results
    tbl_results_soft = Table.read(dir_cstack_output+'/all_soft_results.fits')
    tbl_results_hard = Table.read(dir_cstack_output+'/all_hard_results.fits')

    # calculate the distribution of count rates for individual objects
    # save them to the following directory

    save_all_prob_dists(tbl_results_soft,dir_distribution+'/soft/')
    save_all_prob_dists(tbl_results_hard,dir_distribution+'/hard/')


    # combining the soft and the hard band distributions to get the full band distribution
    # all_band_distribution()

    # get the distribution of stacked photon rates
    get_master_distribution_new('soft')
    get_master_distribution_new('hard')
    # get_master_distribution_new('all')
