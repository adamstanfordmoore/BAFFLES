#!/usr/bin/env python3
"""
Adam Stanford-Moore,Eric Nielsen, Bruce Macintosh, Rob De Rosa, Ian Czekala
Stanford University Physics Department
5/22/20
BAFFLES: Bayesian Ages for Field LowEr-mass Stars
"""

import ca_constants as const
from scipy import interpolate
from scipy.stats import norm,lognorm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import probability as prob
import fitting as my_fits
import numpy as np
import plotting as my_plot
import bisect
import sys
import copy
import utils
import time

# shortcut to quickly computing using default grids the posteriors for calcium and/or lithium
def baffles_age(bv=None,rhk=None,li=None,bv_err=None,li_err = None,upperLim=False,
            maxAge=None,fileName='baffles',pdfPage=None,showPlots=True,
            savePlots=False,savePostAsText=False):
    if (not rhk and not li):
        raise RuntimeError("Must provide at least one of calcium logR'HK or lithium EW")
    if (li and not bv):
        raise RuntimeError("Must provide B-V value with lithium EW")
    
    if (not pdfPage and savePlots):
        pdfPage = PdfPages(fileName + '.pdf')
        
    p = None
    if (rhk):
        baf = age_estimator('calcium')
        p = baf.get_posterior(bv,rhk,pdfPage,showPlots,bv_err,li_err,
                upperLim=None,maxAge=maxAge,mamajekAge=utils.getMamaAge(rhk))
        if (savePostAsText):
            np.savetxt(fileName + "_calcium.csv", list(zip(const.AGE,p.array)), delimiter=",")
        print("Ca Median Age: %.3g Myr, 68%% CI: %.3g - %.3g Myr, 95%% CI: %.3g - %.3g Myr" \
               % (p.stats[2],p.stats[1],p.stats[3],p.stats[0],p.stats[4]))
    p2 = None
    if (li):
        baf2 = age_estimator('lithium')
        p2 = baf2.get_posterior(bv,li,pdfPage,showPlots,bv_err,li_err,
                                upperLim=upperLim,maxAge=maxAge)
        if (savePostAsText):
            np.savetxt(fileName + "_lithium.csv", list(zip(const.AGE,p2.array)),
                       delimiter=",")
        
        if p2.upperLim: print("1 sig lower-lim: %.3g Myr, 2 sig lower-lim: \
            %.3g Myr, 3 sig: %.3g Myr" % (p2.stats[2],p2.stats[1],p2.stats[0]))
        else: print("Li Median Age: %.3g Myr, 68%% CI: %.3g - %.3g Myr, 95%% CI: \
        %.3g - %.3g Myr" % (p2.stats[2],p2.stats[1],p2.stats[3],p2.stats[0],p2.stats[4]))

    p3 = None
    if (p and p2):
        title = 'Final Posterior'
        y = p.array * p2.array
        prob.normalize(const.AGE,y)
        stats = prob.stats(const.AGE,y)
        my_plot.posterior(const.AGE,y,prob.stats(const.AGE,y),title,pdfPage,showPlots)
        p3 = posterior()
        p3.array,p3.stats = y,stats
        print("Final Median Age: %.3g Myr, 68%% CI: %.3g - %.3g, 95%% CI: %.3g - %.3g" \
               % (stats[2],stats[1],stats[3],stats[0],stats[4]))

        if (savePostAsText):
            np.savetxt(fileName + "_product.csv", list(zip(const.AGE,y)), delimiter=",")

    if pdfPage:
        pdfPage.close()

    if p3 is not None: return p3
    return p if p is not None else p2


class posterior:
    def __init__(self):
        self.stats = None  #holds age at different CDF values
        self.array = None  # posterior array
        self.upperLim = False  # if its an upper-limit
        self.stars_posteriors = None #array of individual stellar posteriors if array is a product


class age_estimator:
    #takes in a metal idicator either 'calcium' or 'lithium' denoting which method to use
    # option to input the grid_median as an array or as string referencing saved .npy files
    def __init__(self,metal,grid_median=None,default_grids=True,load_pdf_fit=True):
        self.metal = metal
        self.grid_median = None
        self.const = utils.init_constants(metal)
        if (grid_median):
            self.set_grids(grid_median)
        elif (default_grids):
            self.set_grids(self.const.DEFAULT_MEDIAN_GRID)
        if load_pdf_fit: #allows refresh.py to make without needing these
            self.pdf_fit,self.cdf_fit = my_fits.fit_histogram(metal,fromFile=True)
            #self.pdf_fit,self.cdf_fit = my_fits.fit_student_t(metal,fromFile=True)
    
    def set_grids(self,grid_median):
        if (type(grid_median) == str and type(grid_median) == str):
            if (grid_median[-4:] != '.npy'):
                grid_median += '.npy'
            self.grid_median = np.load(grid_median)
        elif (type(grid_median) == 'numpy.ndarray'):
            self.grid_median = grid_median

    #Takes in bv the (B-V)o corrected color and the metallicity to return a posterior object.
    #Metallicity: log(R'HK) if refering to calcium. log equivalent width per mA if lithium.
    def get_posterior(self,bv,metallicity,pdfPage=None,showPlot=False,\
            bv_uncertainty=None,measure_err = None,upperLim=False,\
            maxAge=None,givenAge=None,givenErr = None,mamajekAge=None,title=None,\
            logPlot = False):
        if bv is None and self.metal=='calcium': bv = 0.65
        assert self.const.BV_RANGE[0] <= bv <= self.const.BV_RANGE[1], \
                "B-V of %.2f out of range. Valid range: " % bv + str(self.const.BV_RANGE)
        if self.metal=='calcium':
            assert self.const.METAL_RANGE[0] <= metallicity <= self.const.METAL_RANGE[1], \
                "Indicator value %.2f out of range. Valid range: " % metallicity + str(self.const.METAL_RANGE)
        elif self.metal=='lithium':
            assert self.const.METAL_RANGE_LIN[0] <= metallicity <= self.const.METAL_RANGE_LIN[1], \
                "Indicator value %.2f out of range. Valid range: " \
                % metallicity + str(self.const.METAL_RANGE_LIN) + " mA"
        
        if mamajekAge == True and self.const.METAL_RANGE[0] <= metallicity <= \
                self.const.METAL_RANGE[1]:
            mamajekAge = utils.getMamaAge(metallicity)

        posterior_arr = self.likelihood(bv,bv_uncertainty,metallicity,measure_err,\
                upperLim) * self.prior(maxAge)
        if all(posterior_arr == 0):
            print("Posterior not well defined. Area is zero so adding constant")
            posterior_arr += 1
        
        prob.normalize(self.const.AGE,posterior_arr)
        
        p_struct = posterior()
        p_struct.array = posterior_arr
        p_struct.stats = prob.stats(self.const.AGE,posterior_arr,upperLim)
        p_struct.upperLim = upperLim
        if (showPlot or pdfPage):
            if (title == None):
                title = '(B-V)o = '+'%.2f' % bv \
                        +', ' + self.metal + ' = %.2f' % metallicity
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,pdfPage,\
                    showPlot,givenAge=givenAge,givenErr=givenErr, mamajekAge=mamajekAge,logPlot=logPlot) 
        return p_struct
    
    def resample_posterior_product(self,bv_arr,metallicity_arr,bv_errs=None,measure_err_arr=None,\
            upperLim_arr=None,maxAge_arr = None, \
            pdfPage=None,showPlot=False,showStars=False,title=None,givenAge=None,givenErr = None,
            sample_num=None,numIter=4):
        if (bv_errs is None):
            bv_errs = [self.const.BV_UNCERTAINTY]*len(bv_arr)
        if upperLim_arr is None:
            upperLim_arr = [False]*len(metallicity_arr)
        if maxAge_arr is None:
            maxAge_arr = [None]*len(metallicity_arr)
        if measure_err_arr is None:
            measure_err_arr = [self.const.MEASURE_ERR]*len(metallicity_arr)
        if sample_num is None: sample_num = 15#len(bv_arr) - 5

        ln_prob = np.zeros(len(self.const.AGE))
        star_post = []
        
        resample_args = (bv_arr,metallicity_arr,bv_errs,measure_err_arr,upperLim_arr,maxAge_arr)
        args = (pdfPage,showPlot,showStars,title,givenAge,givenErr)

        post = prob.resample(self.posterior_product,resample_args,args,sample_num,numIter)

        prob.normalize(self.const.AGE,post)
        p_struct = posterior()
        p_struct.array = post
        p_struct.stats = prob.stats(self.const.AGE,post)

        if (showPlot or pdfPage):
            title = title if title else 'Resampled Posterior Product Age Distribution'
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,\
                    pdfPage,showPlot,star_post,givenAge,givenErr=givenErr)

        return p_struct
    
    def posterior_product(self,bv_arr,metallicity_arr,bv_errs=None,measure_err_arr=None,\
            upperLim_arr=None,maxAge_arr = None, \
            pdfPage=None,showPlot=False,showStars=False,title=None,givenAge=None,givenErr = None):
        if (bv_errs is None):
            bv_errs = [self.const.BV_UNCERTAINTY]*len(bv_arr)
        if upperLim_arr is None:
            upperLim_arr = [False]*len(metallicity_arr)
        if maxAge_arr is None:
            maxAge_arr = [None]*len(metallicity_arr)
        if measure_err_arr is None:
            measure_err_arr = [self.const.MEASURE_ERR]*len(metallicity_arr)
        ln_prob = np.zeros(len(self.const.AGE))
        star_post = []
        
        np.seterr(divide = 'ignore')
        
        if self.metal == 'lithium' and np.mean(metallicity_arr) < 3: 
            metallicity_arr = np.power(10,metallicity_arr) 

        sec_per_star = 0.5
        start=time.time()
        for i in range(len(bv_arr)):
            star_time=time.time()
            y = self.likelihood(bv_arr[i],bv_errs[i],metallicity_arr[i],measure_err_arr[i],\
                  upperLim_arr[i]) * self.prior(maxAge_arr[i])
            prob.normalize(self.const.AGE,y)
            inds = np.nonzero(y)[0]
            ln_prob += np.log(y)
            if (showStars):
                star_post.append(y)
            utils.progress_bar(float(i+1)/len(bv_arr),int((len(bv_arr)-(i+1))*sec_per_star))
            #exp moving average
            sec_per_star = sec_per_star + 0.1*((time.time() - star_time) - sec_per_star)
            #sec_per_star = (time.time() - start)/float(i+1)

        print("Finished %d stars. Average time per star: %.2f seconds." \
              % (len(bv_arr),(time.time() - start)/len(bv_arr)))

        post = np.exp(ln_prob - np.max(ln_prob))  #prevent underflow
        prob.normalize(self.const.AGE,post)
        p_struct = posterior()
        p_struct.array = post
        p_struct.stats = prob.stats(self.const.AGE,post)
        p_struct.stars_posteriors = star_post

        if (showPlot or pdfPage):
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,\
                              pdfPage,showPlot,star_post,givenAge,givenErr=givenErr,\
                              bv_arr = bv_arr,metal=self.metal)
        return p_struct
   
    # Prior on age
    def prior(self,maxAge=None):
        agePrior = 1
        if maxAge is not None and maxAge < self.const.GALAXY_AGE:
            agePrior = self.const.AGE <= maxAge
        return agePrior
    
    def calcium_likelihood(self,bv,rhk):
        mu = self.grid_median 
        #like gaussian likelihood, divide by std which scales height of function
        pdfs = self.pdf_fit(rhk - mu)
        assert (pdfs >= 0).all(), "Error in numerical fit_histogram" + str(pdfs)
        return np.sum(pdfs,axis=0)

    # axes are (axis0=BV,axis1=AGE,axis2=Li)
    # assumes li is linear space
    def likelihood(self,bv,bv_uncertainty,li,measure_err,isUpperLim):
        if self.metal == 'calcium': return self.calcium_likelihood(bv,li)
        if not bv_uncertainty:
            bv_uncertainty = self.const.BV_UNCERTAINTY
        if not measure_err:
            measure_err = self.const.MEASURE_ERR
        
        num_points = self.const.NUM_BV_POINTS if bv_uncertainty <= 0.03 else self.const.NUM_BV_POINTS * (bv_uncertainty/0.03)
        
        # By spacing the points according to area, the weighting for each B-V point is 
        #implicit.  More points are clustered near bv than farther
        BV = prob.gaussian_cdf_space(bv,bv_uncertainty,num_points, sig_lim=3)
        
        f = interpolate.interp2d(self.const.AGE,self.const.BV_S,self.grid_median)
        mu = f(self.const.AGE,BV)

        if isUpperLim:
            #integration done in logspace with log li and log mu
            astro_gauss = self.cdf_fit(np.log10(li) - mu)
            final_sum = np.sum(astro_gauss,axis=0)
            return final_sum

        mu = mu.reshape(mu.shape[0],mu.shape[1],1)

        li_gauss = prob.gaussian(self.const.METAL,li,measure_err)
        mask = li_gauss > prob.FIVE_SIGMAS
        li_gauss = li_gauss[mask]
        li_gauss = li_gauss.reshape(1,1,len(li_gauss))
        METAL = self.const.METAL[mask]
        METAL = METAL.reshape(1,1,len(METAL))
        
        astro_gauss = self.pdf_fit(np.log10(METAL) - mu)/METAL
        product = li_gauss*astro_gauss
        integral = np.trapz(product,METAL,axis=2) # now 2d matrix
        final_sum = np.sum(integral,axis=0)
        return final_sum

    #calculates and returns a 2D array of median b-v and age
    #omit_cluster specifies a cluster index to remove from the fits to make the grids without a cluster
    def make_grids(self,bv_li,fits,upper_lim=None,medianSavefile=None,\
            setAsDefaults=False, omit_cluster=None):
        if (medianSavefile and medianSavefile[-4:] == '.npy'):
            medianSavefile = medianSavefile[:-4]

        primordial_li_fit = None
        bldb_fit = None
        if (self.metal == 'lithium'):
            primordial_li_fit = my_fits.MIST_primordial_li()
            bldb_fit = my_fits.bldb_fit(fits)
        
        median_rhk, sigma = [],[]
        for bv in self.const.BV_S:
            rhk,scatter,CLUSTER_AGES,_ = my_fits.get_valid_metal(bv,fits,self.const,
                    primordial_li_fit,bldb_fit,omit_cluster)
            mu,sig,_ = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,self.metal,
                                           omit_cluster)
            median_rhk.append(mu(self.const.AGE))
            
        median_rhk = np.array(median_rhk)

        if medianSavefile:
            np.save(medianSavefile, median_rhk)

        if (setAsDefaults):
            self.set_default_grids(medianSavefile)
    
        self.grid_median = median_rhk

    #given an x_value, which is the location to evaluate and an array of fits,
    #it calculates the y-value at x_val for each fit and returns an array of them
    def arrayFromFits(self,x_val,fits):
        arr = []
        for fit in fits:
            arr.append(fit(x_val))
        return arr

    def get_grids(self):
        return self.grid_median

    #changes constant file to have these names as default grid names
    def set_default_grids(self,grid_median):
        filename = 'ca_constants.py'
        if (self.metal=='lithium'):
            filename = 'li_constants.py'
        lines = open(filename,'r').readlines()
        for l in range(len(lines)):
            if (lines[l][:14] == 'DEFAULT_MEDIAN'):
                lines[l] = 'DEFAULT_MEDIAN_GRID = "' + grid_median + '.npy"\n'
        out = open(filename,'w')
        out.writelines(lines)
        out.close()

if  __name__ == "__main__":
    const = utils.init_constants('lithium')
    err = "Usage:  python baffles.py -bmv <B-V> -rhk <Log10(R\'HK)> -li <EWLi> [-bmv_err <> -li_err <> -ul -maxAge <13000> -plot -s -savePlot -filename <> -help]"
    
    help_msg = "\n\
    -bmv corrected (B-V)o of the star (optional for calcium)\n\
    -rhk <> log base 10 of the R'HK value \n\
    -li <> EW measure in milli-angstroms (=0.1 pm) \n\
    \noptional flags:\n\
    -ul indicates that the log(EW/mA) reading is an upper-limit reading. \n\
    -maxAge allows user to input max posible age of star (Myr) if upper-limit flag is used. default is %d \n\
    -bmv_err <float uncertainity> provides the uncertainty on B-V with default %.2f \n\
    -li_err <float uncertainity> provides the uncertainty in LiEW measurement with default %dmA \n\
    -s saves posteriors as .csv as age, probability in two 1000 element columns.\n\
    -plot plots and shows the PDF. \n\
    -savePlot saves the plotted posteriors to a pdf. \n\
    -filename <name of file w/o extension> name of files to be saved: name.pdf is graphs, name_calcium.csv/name_lithium.csv/name_product.csv are posterior csv files for calcium/lithium/product respectively.  \n\
    -help prints this message \n" % (const.GALAXY_AGE,const.BV_UNCERTAINTY, const.MEASURE_ERR)
    
    argv = sys.argv #starts with 'baffles.py'
    if (len(argv) < 3 or '-help' in argv):
        print(err,help_msg)
        sys.exit()
    bv = None
    bv_err,li_err = None,None
    rhk, li = None,None
    save = False
    showPlots = False
    fileName = 'baffles'
    savePlots = False
    upperLim = False
    maxAge = const.GALAXY_AGE
    valid_flags = ['-bmv','-rhk','-li','-li_err','-bmv_err','-plot','-savePlot','-ul','-maxAge','-s','-filename','-help']
    extra_flags = ['-Plot','-plots','-Plots','-savePlots','-saveplots','-saveplot','-UL','-save']
    for i,ar in enumerate(argv[1:]):
        if ar not in valid_flags and ar not in extra_flags \
            and not utils.isFloat(ar) and argv[i] != '-filename':
            print("Invalid flag '" + ar + "'. Did you mean one of these:")
            print(valid_flags)
            exit()
    try:
        if ('-bmv' in argv):
            bv = float(argv[argv.index('-bmv') + 1])
        if ('-rhk' in argv):
            rhk = float(argv[argv.index('-rhk') + 1])
            import ca_constants as const
            if bv is not None and (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print("B-V out of range. Must be in range " + str(const.BV_RANGE))
                sys.exit()
            if (not (const.METAL_RANGE[0] <= rhk <= const.METAL_RANGE[1])):
                print("Log(R\'HK) out of range. Must be in range " + str(const.METAL_RANGE))
                sys.exit()
        if ('-li' in argv):
            li = float(argv[argv.index('-li') + 1])
            import li_constants as const
            if 0 < li < 3:
                print("Interpretting LiEW as log(LiEW)")
                li = 10**li
            
            if (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print("B-V out of range. Must be in range " + str(const.BV_RANGE))
                sys.exit()
            if (not (const.METAL_RANGE_LIN[0] <= li <= const.METAL_RANGE_LIN[1])):
                print("Li EW out of range. Must be in range " + str(const.METAL_RANGE) + " mA")
                sys.exit()
        if ('-s' in argv or '-save' in argv):
            save = True
        if ('-plot' in argv or '-Plot' in argv or '-plots' in argv or '-Plots' in argv):
            showPlots = True
        if ('-saveplot' in argv or '-saveplots' in argv or '-savePlots' in argv or '-savePlot' in argv):
            savePlots = True
        if ('-filename' in argv):
            fileName = argv[argv.index('-filename') + 1]
        if ('-bmv_err' in argv):
            bv_err = float(argv[argv.index('-bmv_err') + 1])
        if ('-li_err' in argv):
            li_err = float(argv[argv.index('-li_err') + 1])
        if ('-ul' in argv or '-UL' in argv):
            upperLim = True
        if ('-maxAge' in argv):
            maxAge = float(argv[argv.index('-maxAge') + 1])
    except IndexError:
        print(err)
    except ValueError:
        print(err)
    
    baffles_age(bv,rhk,li,bv_err,li_err,upperLim,maxAge,fileName,showPlots=showPlots,
                savePlots=savePlots, savePostAsText=save)
