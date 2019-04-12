#!/usr/bin/env python2.7
"""
Adam Stanford-Moore,Eric Nielsen, Bruce Macintosh, Rob De Rosa
Stanford University Physics Department
8/28/18
BAFFLES: Bayesian Ages for Field LowEr-mass Stars
"""

import ca_constants as const
BV_UNCERTAINTY = .002
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import probability as prob
import fitting as my_fits
import numpy as np
import plotting as my_plot
import bisect
import sys
import copy
# shortcut to quickly computing using default grids the posteriors for calcium and/or lithium
def baffles_age(bv,rhk=None,li=None,bv_err=BV_UNCERTAINTY,fileName='baffles',pdfPage=None,showPlots=True,noPlots=False,savePostAsText=False):
    if (not rhk and not li):
        raise RuntimeError("Must provide at least one of calcium logR'HK or lithium log(EW/mA)")
    
    if (not pdfPage and not noPlots):
        pdfPage = PdfPages(fileName + '.pdf')
        
    p = None
    if (rhk):
        baf = age_estimator('calcium')
        p = baf.get_posterior(bv,rhk,pdfPage,showPlots,bv_uncertainty=bv_err,mamajekAge=getMamaAge(rhk))
        if (savePostAsText):
            np.savetxt(fileName + "_calcium.csv", zip(const.AGE,p.array), delimiter=",")
    p2 = None
    if (li):
        baf2 = age_estimator('lithium')
        p2 = baf2.get_posterior(bv,li,pdfPage,showPlots,bv_uncertainty=bv_err)
        if (savePostAsText):
            np.savetxt(fileName + "_lithium.csv", zip(const.AGE,p2.array), delimiter=",")

    if (p and p2):
        title = ' Calcium/Lithium Posterior Product'
        y = p.array * p2.array
        prob.normalize(const.AGE,y)
        my_plot.posterior(const.AGE,y,prob.stats(const.AGE,y),title,pdfPage,showPlots)

        if (savePostAsText):
            np.savetxt(fileName + "_product.csv", zip(const.AGE,y), delimiter=",")


class posterior:
    def __init__(self):
        self.stats = None
        self.array = None


class age_estimator:
    #takes in a metal idicator either 'calcium' or 'lithium' denoting which method to use
    # option to input the grid_median and grid_sigma as arrays or as strings referencing saved .npy files
    def __init__(self,metal,grid_median=None,grid_sigma=None,default_grids=True):
        self.metal = metal
        self.grid_median = None
        self.grid_sigma = None
        self.const = self.init_constants(metal)
        if (grid_median and grid_sigma):
            self.set_grids(grid_median,grid_sigma)
        elif (default_grids):
            self.set_grids(self.const.DEFAULT_MEDIAN_GRID,self.const.DEFAULT_SIGMA_GRID)
    
    def set_grids(self,grid_median,grid_sigma):
        if (type(grid_median) == str and type(grid_median) == str):
            if (grid_median[-4:] != '.npy'):
                grid_median += '.npy'
            if (grid_sigma[-4:] != '.npy'):
                grid_sigma += '.npy'
            self.grid_median, self.grid_sigma = np.load(grid_median), np.load(grid_sigma)
        elif (type(grid_median) == 'numpy.ndarray'):
            self.grid_median, self.grid_sigma = grid_median, grid_sigma

    #Takes in bv the (B-V)o corrected color and the metallicity to return a posterior object.
    #Metallicity: log(R'HK) if refering to calcium. log equivalent width per mA if lithium.
    def get_posterior(self,bv,metallicity,pdfPage=None,showPlot=False,givenAge=None,givenErr = None,bv_uncertainty=BV_UNCERTAINTY,mamajekAge=None,title=None):
        array = self.age_dist_uncertainty(bv,bv_uncertainty,metallicity)
        prob.normalize(self.const.AGE,array)
        p_struct = posterior()
        p_struct.array = array
        p_struct.stats = prob.stats(self.const.AGE,array)

        if (showPlot or pdfPage):
            if (title == None):
                title = 'Posterior Age Probabilty for (B-V)o = '+'%.2f' % bv +', ' + self.metal + ' = %.2f' % metallicity
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,pdfPage,showPlot,givenAge=givenAge,givenErr=givenErr, mamajekAge=mamajekAge) 
        return p_struct
    
    def posterior_product(self,bv_arr,metallicity_arr,pdfPage=None,showPlot=False,showStars=False,title=None,givenAge=None,givenErr = None,bv_errs=None):
        if (not bv_errs):
            bv_errs = [BV_UNCERTAINTY]*len(bv_arr)
        ln_prob = np.zeros(len(self.const.AGE))
        star_post = []
        for i in range(len(bv_arr)):
                y = self.age_dist_uncertainty(bv_arr[i],bv_errs[i],metallicity_arr[i])
                ln_prob += np.log10(y)
                if (showStars):
                    prob.normalize(self.const.AGE,y)
                    star_post.append(y)
        post = np.exp(ln_prob)
        prob.normalize(self.const.AGE,post)
        p_struct = posterior()
        p_struct.array = post
        p_struct.stats = prob.stats(self.const.AGE,post)

        if (showPlot or pdfPage):
                title = title if title else 'Posterior Product Age Distribution'
                my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,pdfPage,showPlot,star_post,givenAge,givenErr=givenErr)

        return p_struct
    
    def init_constants(self,metal):
        if (metal[0].lower() == 'c'):
            self.metal = 'calcium'
            import ca_constants as const
        elif (metal[0].lower() == 'l'):
            self.metal = 'lithium'
            import li_constants as const
        else:
            raise RuntimeError("No metal specified. Please enter lithium or calcium")
        return const

    #returns an array representing the posterior
    def age_dist_uncertainty(self,bv,sig,rhk):
        gauss = prob.gaussian(bv,sig,self.const.BV)
        arr = np.zeros(len(self.grid_median[0]))
        for i in range(len(self.const.BV)):
            if (gauss[i] > prob.FIVE_SIGMAS):
                arr += gauss[i]*self.likelihood(i,rhk)
        return arr

    def age_dist(self,bv,rhk): #returns an array
        b = bisect.bisect_left(self.const.BV,bv)
        return self.likelihood(b,rhk)

    #calculated likelihood given b, the index of BV for the calculation
    def likelihood(self,b,rhk):
        mu = np.array(self.grid_median[b])
        sig = np.array(self.grid_sigma[b])
        return prob.gaussian(mu,sig,rhk)

    #calculates and returns a 2D array of sigma b-v and age
    #omit_cluster specifies a cluster index to remove from the fits to make the grids without a cluster
    def make_grids(self,bv_li,fits,medianSavefile=None,sigmaSavefile=None,setAsDefaults=False, omit_cluster=None,upper_lim=None):
        if (medianSavefile and medianSavefile[-4:] == '.npy'):
            medianSavefile = medianSavefile[:-4]
        if (sigmaSavefile and sigmaSavefile[-4:] == '.npy'):
            sigmaSavefile = sigmaSavefile[:-4]

        primordial_li_fit = None
        ldb_fit = None
        li_scatter_fit = None #fit as function of LiEW
        ca_scatter = None
        if (self.metal == 'lithium'):
            #ngc2264_fit = fits[const.CLUSTER_NAMES.index("NGC2264")][0]
            primordial_li_fit = my_fits.primordial_li()
            bldb_fit = my_fits.bldb_fit(fits)
            li_scatter_fit = my_fits.fit_two_scatters(bv_li,fits,upper_lim=upper_lim,omit_cluster=omit_cluster)
        elif (self.metal == 'calcium'):
            ca_scatter = my_fits.total_scatter(bv_li,fits,omit_cluster)
        
        median_rhk, sigma = [],[]
        for bv in self.const.BV:
            rhk,scatter,CLUSTER_AGES = [],[],[]
            for i in range(len(fits)):
                if (omit_cluster and i == omit_cluster):
                    continue
                r = fits[i][0](bv)
                if (self.const.METAL_RANGE[0] <  r < self.const.METAL_RANGE[1]):
                    rhk.append(r)
                    scatter.append(fits[i][1](bv))
                    CLUSTER_AGES.append(self.const.CLUSTER_AGES[i])
                    #plt.scatter(self.const.CLUSTER_AGES[i],r,label=self.const.CLUSTER_NAMES[i])
            #plt.legend()        
            #plt.show()
            #info added from depletion boundary and primordial lithium
            if (self.metal == 'lithium'):
                CLUSTER_AGES.append(bldb_fit(bv)) 
                rhk.append(self.const.ZERO_LI)
                CLUSTER_AGES.append(self.const.PRIMORDIAL_LI_AGE)
                rhk.append(primordial_li_fit(bv))

            fit = None
            if (self.metal == 'calcium'):
                fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,2)
            elif (self.metal == 'lithium'):
                if (.76 <= bv <= .94): #Patch to fix bad region. NEED TO FIX BETTER
                    fit = my_fits.piecewise([0,2.2,3],[2.5,2.1,.5])
                else:
                    fit = my_fits.constrained_poly_fit(np.log10(CLUSTER_AGES),rhk,0)
            """
            if (len(rhk) == 1):
                f = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,0)
            if (len(rhk) == 2):
                f = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,1)
            elif (self.metal == 'calcium'):
                f = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,2) #like mamajek polynomial is non-linear
            elif (self.metal == 'lithium'):   
                if (bv <=1.35):
                    f = interpolate.interp1d(np.log10(CLUSTER_AGES),rhk, fill_value='extrapolate')
                elif (bv <= 1.55):
                    f = my_fits.constrained_poly_fit(np.log10(CLUSTER_AGES),rhk,0)
                    #plt.semilogx(CLUSTER_AGES,rhk)
                    #plt.semilogx(self.const.AGE,f(np.log10(self.const.AGE)))
                    #plt.legend()        
                    #plt.show()
                else:
                    f = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,1)
            """
            median_rhk.append(fit(np.log10(self.const.AGE))) #uses determined function 

            #4 different methods for handling scatter -- total detrended mean, mean clusters, best-fit, linear interp
            if (self.metal == 'calcium'):
                sigma.append([ca_scatter for i in range(len(self.const.AGE))])
            elif (self.metal == 'lithium'):
                sig = li_scatter_fit(fit(np.log10(self.const.AGE)))
                sigma.append(sig)
            
            #m = np.mean(scatter)
            #sigma.append([m for i in range(len(self.const.AGE))])
            #sigma.append(np.poly1d(np.polyfit(np.log10(self.const.CLUSTER_AGES), scatter, 1))(np.log10(self.const.AGE)))
            #g = interpolate.interp1d(np.log10(self.const.CLUSTER_AGES),scatter, fill_value='extrapolate')
            #sigma.append(g(np.log10(self.const.AGE))) #linear extrapolate

        if (medianSavefile and sigmaSavefile):
            np.save(medianSavefile, median_rhk)
            np.save(sigmaSavefile, sigma)

        if (setAsDefaults):
            self.set_default_grids(medianSavefile,sigmaSavefile)
    
        self.grid_median, self.grid_sigma = median_rhk, sigma

    #given an x_value, which is the location to evaluate and an array of fits,
    #it calculates the y-value at x_val for each fit and returns an array of them
    def arrayFromFits(self,x_val,fits):
        arr = []
        for fit in fits:
            arr.append(fit(x_val))
        return arr

    def get_grids(self):
        return self.grid_median, self.grid_sigma

    #changes constant file to have these names as default grid names
    def set_default_grids(self,grid_median,grid_sigma):
        filename = 'ca_constants.py'
        if (self.metal=='lithium'):
            filename = 'li_constants.py'
        lines = open(filename,'r').readlines()
        for l in range(len(lines)):
            if (lines[l][:14] == 'DEFAULT_MEDIAN'):
                lines[l] = 'DEFAULT_MEDIAN_GRID = "' + grid_median + '.npy"\n'
            if (lines[l][:13] == 'DEFAULT_SIGMA'):
                lines[l] = 'DEFAULT_SIGMA_GRID = "' + grid_sigma + '.npy"\n'
        out = open(filename,'w')
        out.writelines(lines)
        out.close()

def getMamaAge(r):
        return np.power(10,-38.053 - 17.912*r - 1.6675*r*r)/1e6

def getMamaRHK(age):
    log_age = np.log10(np.array(age)*1e6)
    return 8.94 - 4.849*log_age + .624*log_age**2 - .028 * log_age**3


if  __name__ == "__main__":
    err = "Usage:  python baffles.py -bv <B-V> -rhk <Log(R\'HK)> -li <Log(EW/mA)> [-bverr <> -noplot -s -showplots -filename <> -help]"
    
    help_msg = "\noptional flags:\n\
    -bverr <float uncertainity> provides the uncertainty on B-V with default .002 \n\
    -noplot suppresses all plotting and just generate .csv files (-s extension is implied. Used with -showplots it prevents saving)\n\
    -s saves posteriors as .csv \n\
    -showplots to show posterior plots for user before saving\n\
    -filename <name of file w/o extension> name of files to be saved: name.pdf is graphs, name_calcium.csv/name_lithium.csv/name_product.csv are posterior csv files for calcium/lithium/product respectively.  csv is stored as age, probability in two 13000 element columns.\n\
    -help prints this message"
    
    argv = sys.argv
    if (len(argv) < 5 or '-help' in argv):
        print err,help_msg
        sys.exit()
    bv = None
    bv_err = BV_UNCERTAINTY
    rhk, li = None,None
    save = False
    showPlots = False
    fileName = 'baffles'
    noPlots = False
    try:
        bv = float(argv[argv.index('-bv') + 1])
        if ('-rhk' in argv):
            rhk = float(argv[argv.index('-rhk') + 1])
            import ca_constants as const
            if (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print "B-V out of range. Must be in range " + str(const.BV_RANGE)
                sys.exit()
            if (not (const.METAL_RANGE[0] <= rhk <= const.METAL_RANGE[1])):
                print "Log(R\'HK) out of range. Must be in range " + str(const.METAL_RANGE)
                sys.exit()
        if ('-li' in argv):
            li = float(argv[argv.index('-li') + 1])
            import li_constants as const
            if (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print "B-V out of range. Must be in range " + str(const.BV_RANGE)
                sys.exit()
            if (not (const.METAL_RANGE[0] <= li <= const.METAL_RANGE[1])):
                print "Log(EW/mA) out of range. Must be in range " + str(const.METAL_RANGE)
                sys.exit()
        if ('-s' in argv):
            save = True
        if ('-showplots' in argv):
            showPlots = True
        if ('-noplot' in argv):
            noPlots = True
            save = True
        if ('-filename' in argv):
            fileName = argv[argv.index('-filename') + 1]
        if ('-bverr' in argv):
            bv_err = float(argv[argv.index('-bverr') + 1])
            
    except IndexError:
        print err
    except ValueError:
        print err
    
    baffles_age(bv,rhk,li,bv_err,fileName,showPlots=showPlots,noPlots=noPlots, savePostAsText=save)
