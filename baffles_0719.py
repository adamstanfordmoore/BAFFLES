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
def baffles_age(bv,rhk=None,li=None,bv_err=BV_UNCERTAINTY,upperLim=False,maxAge=const.GALAXY_AGE,fileName='baffles',pdfPage=None,showPlots=True,noPlots=False,savePostAsText=False):
    if (not rhk and not li):
        raise RuntimeError("Must provide at least one of calcium logR'HK or lithium log(EW/mA)")
    
    if (not pdfPage and not noPlots):
        pdfPage = PdfPages(fileName + '.pdf')
        
    p = None
    if (rhk):
        baf = age_estimator('calcium')
        p = baf.get_posterior(bv,rhk,pdfPage,showPlots,bv_uncertainty=bv_err,isUpperLim=upperLim,maxAge=maxAge,mamajekAge=my_fits.getMamaAge(rhk))
        if (savePostAsText):
            np.savetxt(fileName + "_calcium.csv", zip(const.AGE,p.array), delimiter=",")
    p2 = None
    if (li):
        baf2 = age_estimator('lithium')
        p2 = baf2.get_posterior(bv,li,pdfPage,showPlots,bv_uncertainty=bv_err,isUpperLim=upperLim,maxAge=maxAge)
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
        self.upperLim = False


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
    def get_posterior(self,bv,metallicity,pdfPage=None,showPlot=False,\
            bv_uncertainty=BV_UNCERTAINTY,measure_err = None,upperLim=False,\
            maxAge=None,givenAge=None,givenErr = None,mamajekAge=None,title=None,\
            logPlot = False):
        assert self.const.BV_RANGE[0] <= bv <= self.const.BV_RANGE[1], \
                "B-V of %.2f out of range. Valid range: " % bv + str(self.const.BV_RANGE)
        assert self.const.METAL_RANGE[0] <= metallicity <= self.const.METAL_RANGE[1], \
                "Indicator value %.2f out of range. Valid range: " % metallicity + str(self.const.METAL_RANGE)
        
        if self.metal == 'lithium': metallicity = 10**metallicity 
        measure_err = 15
        
        posterior_arr = self.likelihood(bv,bv_uncertainty,metallicity,measure_err,\
                upperLim) * self.prior(maxAge)
        prob.normalize(self.const.AGE,posterior_arr)
        
        p_struct = posterior()
        p_struct.array = posterior_arr
        p_struct.stats = prob.stats(self.const.AGE,posterior_arr,upperLim)
        p_struct.upperLim = upperLim
        if (showPlot or pdfPage):
            if (title == None):
                title = 'Posterior Age Probabilty for (B-V)o = '+'%.2f' % bv \
                        +', ' + self.metal + ' = %.2f' % metallicity
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,pdfPage,\
                    showPlot,givenAge=givenAge,givenErr=givenErr, mamajekAge=mamajekAge,logPlot=logPlot) 
        return p_struct
    
    def resample_posterior_product(self,bv_arr,metallicity_arr,bv_errs=None,measure_err_arr=None,\
            upperLim_arr=None,maxAge_arr = None, \
            pdfPage=None,showPlot=False,showStars=False,title=None,givenAge=None,givenErr = None,
            sample_num=None,numIter=4):
        if (bv_errs is None):
            bv_errs = [BV_UNCERTAINTY]*len(bv_arr)
        if upperLim_arr is None:
            upperLim_arr = [False]*len(metallicity_arr)
        if maxAge_arr is None:
            maxAge_arr = [None]*len(metallicity_arr)
        if measure_err_arr is None:
            measure_err_arr = [15]*len(metallicity_arr)
        if sample_num is None: sample_num = 15#len(bv_arr) - 5

        ln_prob = np.zeros(len(self.const.AGE))
        star_post = []
        
        measure_err = 15
        
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
            bv_errs = [BV_UNCERTAINTY]*len(bv_arr)
        if upperLim_arr is None:
            upperLim_arr = [False]*len(metallicity_arr)
        if maxAge_arr is None:
            maxAge_arr = [None]*len(metallicity_arr)
        if measure_err_arr is None:
            measure_err_arr = [15]*len(metallicity_arr)
        ln_prob = np.zeros(len(self.const.AGE))
        star_post = []
        
        
        measure_err = 15
        if self.metal == 'lithium': metallicity_arr = np.power(10,metallicity_arr) 
        
        sec_per_star = 1
        start=time.time()
        for i in range(len(bv_arr)):
            star_time=time.time()
            y = self.likelihood(bv_arr[i],bv_errs[i],metallicity_arr[i],measure_err_arr[i],\
                  upperLim_arr[i]) * self.prior(maxAge_arr[i])
            prob.normalize(self.const.AGE,y)
            #print "Age: %.3g / mamaAge: %.3g" % (prob.stats(self.const.AGE,y,upperLim_arr[i])[2], utils.getMamaAge(metallicity_arr[i]))
            #y = self.age_dist_uncertainty(bv_arr[i],bv_errs[i],metallicity_arr[i],upperLim_arr[i],\
            #    maxAge_arr[i])
            #if np.any(y <= 0):
            #    print "\n",bv_arr[i],metallicity_arr[i]
            #    print self.const.AGE[y <= 0][0]
            #    #plt.plot(self.const.AGE,y)
            #    #plt.show()
            ln_prob += np.log(y)
            if (showStars):
                #prob.normalize(self.const.AGE,y)
                star_post.append(y)
            utils.progress_bar(float(i+1)/len(bv_arr),int((len(bv_arr)-(i+1))*sec_per_star))
            sec_per_star = sec_per_star + 0.2*((time.time() - star_time) - sec_per_star)
            #sec_per_star = (time.time() - start)/float(i+1)

        print "Finished %d stars. Average time per star: %.2f seconds." % (len(bv_arr),(time.time() - start)/len(bv_arr))

        post = np.exp(ln_prob - np.max(ln_prob))
        prob.normalize(self.const.AGE,post)
        p_struct = posterior()
        p_struct.array = post
        p_struct.stats = prob.stats(self.const.AGE,post)

        if (showPlot or pdfPage):
            title = title if title else 'Posterior Product Age Distribution'
            my_plot.posterior(self.const.AGE, p_struct.array, p_struct.stats,title,\
                    pdfPage,showPlot,star_post,givenAge,givenErr=givenErr,bv_arr = bv_arr,metal=self.metal)

        #return p_struct.array
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
    def age_dist_uncertainty(self,bv,sig,rhk,isUpperLim=False,maxAge=None):
        gauss = norm.logpdf(self.const.BV,loc=bv,scale=sig)
        arr = np.zeros(len(self.grid_median[0]))
        a = -5
        for i in range(len(self.const.BV)):
            if (gauss[i] > np.log(prob.FIVE_SIGMAS)):
                arr += np.exp(gauss[i] + self.log_likelihood(i,rhk,isUpperLim)-a)
        agePrior = 1 if not maxAge else self.const.AGE <= maxAge
        return arr*np.exp(a)*agePrior

    #def age_dist(self,bv,rhk): #returns an array
    #    b = bisect.bisect_left(self.const.BV,bv)
    #    return self.likelihood(b,rhk)

    #calculated likelihood given b, the index of BV for the calculation
    def log_likelihood(self,b,rhk,isUpperLim):
        mu = np.array(self.grid_median[b])
        sig = np.array(self.grid_sigma[b])
        assert not utils.hasNan(mu) and not utils.hasNan(sig),"Found Nan in fit!"
        if isUpperLim:
            return norm.logcdf(rhk,loc=mu,scale=sig)
        return prob.log_gaussian(rhk,mu,sig)
    
    def calcium_likelihood(self,bv,rhk):
        b = 0 #bisect.bisect_left(self.const.BV,bv) ######### SHORTCUT ##########
        mu = np.array(self.grid_median[b])
        sig = np.array(self.grid_sigma[b])
        return prob.gaussian(rhk,mu,sig)


    # Prior on age
    def prior(self,maxAge=None):
        agePrior = 1 if maxAge is None else self.const.AGE <= maxAge
        return agePrior
    
    # axes are (axis0=BV,axis1=AGE,axis2=Li)
    def likelihood(self,bv,bv_uncertainty,li,measure_err,isUpperLim):
        if self.metal == 'calcium': return self.calcium_likelihood(bv,li)
        start = time.time()

        pdf_fit,cdf_fit = my_fits.fit_histogram(metal=self.const.METAL_NAME,fromFile=True)
        #pdf_fit = lambda x: norm.pdf(x,loc=0,scale=0.17)
        #cdf_fit = lambda x: norm.cdf(x,loc=0,scale=0.17)
        #sig = 0.17 #self.astrophysical_scatter#np.array(self.grid_sigma[b])
        
        #b = bisect.bisect_left(self.const.BV,bv)
        BV = np.linspace(np.max(bv - 4*bv_uncertainty,self.const.BV[0]),\
                np.min(bv + 4*bv_uncertainty,self.const.BV[1]),15)
        bv_gauss = prob.gaussian(BV,bv,bv_uncertainty)
        #mask = bisect.bisect_left(self.const.BV,bv)#bv_gauss > prob.FIVE_SIGMAS
        #bv_gauss = bv_gauss[mask]
        #bv_gauss = bv_gauss.reshape(1,1,1)
        
        #plt.plot(BV,bv_gauss)
        #plt.show()
        bv_gauss = bv_gauss.reshape(len(bv_gauss),1,1)
        #print "BV-len",bv_gauss.shape
        #BV = self.const.BV[mask]
        #BV = BV.reshape(len(BV),1,1)

        #print "Median shape",self.grid_median.shape
        #mu = np.array(self.grid_median[mask]) #now of B-V and age
        #print self.grid_median.shape
        f = interpolate.interp2d(self.const.AGE,self.const.BV_S,self.grid_median)
        mu = f(self.const.AGE,BV)

        if isUpperLim:
            #integration done in logspace with log li and log mu
            #rv = norm(loc=mu,scale=sig)
            #p_arr = np.sum(rv.cdf(np.log10(li)),axis=0)
            #return prob.normalize(self.const.AGE,p_arr)
            astro_gauss = cdf_fit(np.log10(li) - mu)
            final_sum = np.sum(astro_gauss,axis=0)
            prob.normalize(self.const.AGE,final_sum)
            #print "UL time %.2f" % (time.time() - start)
            return final_sum

        #mu = mu.reshape(1,len(mu),1)  #CHANGED 10**
        mu = mu.reshape(mu.shape[0],mu.shape[1],1)  #CHANGED 10**

        #start_sm = time.time()
        li_gauss = prob.gaussian(self.const.METAL,li,measure_err)
        #print "li-gauss shape",li_gauss.shape
        mask = li_gauss > prob.FIVE_SIGMAS
        li_gauss = li_gauss[mask]
        li_gauss = li_gauss.reshape(1,1,len(li_gauss))
        #print "new shape",li_gauss.shape
        METAL = self.const.METAL[mask]
        METAL = METAL.reshape(1,1,len(METAL))
        #print METAL.shape
        #print "mu",mu.shape

        #start = time.time()
        #shifted_metal = METAL - mu
        #print shifted_metal.shape
        astro_gauss = pdf_fit(np.log10(METAL) - mu)/METAL
        #print astro_gauss.shape,'time',time.time() - start
        #print astro_gauss.shape

        #y = METAL/10**mu
        #astro_gauss2 = prob.lognorm(y,s=0.17) / 10**mu

        #print 'Diff', np.sum(np.abs(astro_gauss-astro_gauss2))/astro_gauss.shape[0]/astro_gauss.shape[1]/astro_gauss.shape[2]
        
        #rv = lognorm(s=sig,scale=mu)
        #astro_gauss = rv.pdf(METAL)
        #print("Astro Time",time.time() - start)
        #start = time.time()
        product = li_gauss*astro_gauss
        #print("Product Time",time.time() - start)
        #print astro_gauss.shape
        
        #start = time.time()
        integral = np.trapz(product,METAL,axis=2)
        #print "integral",integral.shape
        final_sum = np.sum(integral,axis=0)
        #prob.normalize(self.const.AGE,final_sum)
        #print final_sum.shape
        #end = time.time()
        #print("Integral Time",end - start)
        #stop_sm = time.time()
        #time_small = stop_sm - start_sm
        #print "Time:", time_small

        #start_no = time.time()
        #sum_noBV = self.get_age_posterior_noBV_err(bv,bv_uncertainty,li,measure_err,isUpperLim)
        #time_noBV = time.time() - start_no

        #plt.plot(METAL[0,0,:],li_gauss[0,0,:],label="measurement gaussian")
        #plt.legend()
        #plt.xlabel("EW mA")
        #plt.show()
        #plt.loglog(self.const.AGE,np.abs(sum_noBV - final_sum)/final_sum,label="Percent Error")
        #plt.title("Time saved:%.3f (%.3f,%.3f)" % (time_small - time_noBV,time_small,time_noBV))
        #plt.ylabel("Percent Error")
        #plt.xlabel("Age")
        #plt.show()

        #plt.semilogx(self.const.AGE,sum_noBV,label="no BV uncertainty")
        #plt.semilogx(self.const.AGE,final_sum,label="BV uncertainty")
        #plt.legend()
        #plt.xlabel("Age")
        #plt.show()
        
        #print "time %.2f" % (time.time() - start)
        return final_sum

    """
    # axes are (axis0=BV,axis1=AGE,axis2=Li)
    def likelihood(self,bv,bv_uncertainty,li,measure_err,isUpperLim):
        if self.metal == 'calcium': return self.calcium_likelihood(bv,li)
        
        sig = 0.17 #self.astrophysical_scatter#np.array(self.grid_sigma[b])
        
        b = bisect.bisect_left(self.const.BV,bv)
        bv_gauss = prob.gaussian(self.const.BV,bv,bv_uncertainty)
        mask = bv_gauss > prob.FIVE_SIGMAS
        bv_gauss = bv_gauss[mask]
        bv_gauss = bv_gauss.reshape(len(bv_gauss),1,1)
        #print "BV-len",bv_gauss.shape
        #BV = self.const.BV[mask]
        #BV = BV.reshape(len(BV),1,1)

        #print "Median shape",self.grid_median.shape
        mu = np.array(self.grid_median[mask]) #now of B-V and age
        if isUpperLim:
            #integration done in logspace with log li and log mu
            rv = norm(loc=mu,scale=sig)
            p_arr = np.sum(rv.cdf(np.log10(li)),axis=0)
            return prob.normalize(self.const.AGE,p_arr)
        
        mu = 10**mu.reshape(mu.shape[0],mu.shape[1],1)
        #print "mu",mu.shape
        
        
        start_sm = time.time()
        li_gauss = prob.gaussian(self.const.METAL,li,measure_err)
        #print "li-gauss shape",li_gauss.shape
        mask = li_gauss > prob.FIVE_SIGMAS
        li_gauss = li_gauss[mask]
        li_gauss = li_gauss.reshape(1,1,len(li_gauss))
        #print "new shape",li_gauss.shape
        METAL = self.const.METAL[mask]
        METAL = METAL.reshape(1,1,len(METAL))
        
        start = time.time()
        y = METAL/mu
        astro_gauss = prob.lognorm(y,s=sig) / mu
        
        #rv = lognorm(s=sig,scale=mu)
        #astro_gauss = rv.pdf(METAL)
        #print("Astro Time",time.time() - start)
        start = time.time()
        product = li_gauss*astro_gauss
        #print("Product Time",time.time() - start)
        #print astro_gauss.shape
        
        start = time.time()
        integral = np.trapz(product,METAL,axis=2)
        #print "integral",integral.shape
        final_sum = np.sum(integral,axis=0)
        #prob.normalize(self.const.AGE,final_sum)
        #print final_sum.shape
        end = time.time()
        #print("Integral Time",end - start)
        stop_sm = time.time()
        time_small = stop_sm - start_sm
        #print "Time:", time_small

        #start_no = time.time()
        #sum_noBV = self.get_age_posterior_noBV_err(bv,bv_uncertainty,li,measure_err,isUpperLim)
        #time_noBV = time.time() - start_no

        #plt.plot(METAL[0,0,:],li_gauss[0,0,:],label="measurement gaussian")
        #plt.legend()
        #plt.xlabel("EW mA")
        #plt.show()
        #plt.loglog(self.const.AGE,np.abs(sum_noBV - final_sum)/final_sum,label="Percent Error")
        #plt.title("Time saved:%.3f (%.3f,%.3f)" % (time_small - time_noBV,time_small,time_noBV))
        #plt.ylabel("Percent Error")
        #plt.xlabel("Age")
        #plt.show()

        #plt.semilogx(self.const.AGE,sum_noBV,label="no BV uncertainty")
        #plt.semilogx(self.const.AGE,final_sum,label="BV uncertainty")
        #plt.legend()
        #plt.xlabel("Age")
        #plt.show()
        
        
        return final_sum
    """

    def get_age_posterior_noBV_err(self,bv,bv_uncertainty,li,measure_err,isUpperLim):
        b = bisect.bisect_left(self.const.BV,bv)
        mu = 10**np.array(self.grid_median[b])
        sig = 0.17 #self.astrophysical_scatter#np.array(self.grid_sigma[b])
        
        if isUpperLim:
            return prob.normalize(self.const.AGE,norm.cdf(li,loc=mu,scale=sig))
        
        #start_full = time.time()
        #li_gauss = prob.gaussian(self.const.METAL,li,measure_err)
        #METAL = self.const.METAL

        #start = time.time()
        #rv = lognorm(s=sig,scale=mu)
        #astro_gauss = rv.pdf(METAL)
        #product = li_gauss*astro_gauss
        #end = time.time()
        #print("Product Time",end - start)
        #print astro_gauss.shape
        
        #start = time.time()
        #integral_full = prob.normalize(self.const.AGE,np.trapz(product,METAL,axis=0))
        #end = time.time()
        #print("Integral Time",end - start)
        #stop_full = time.time()
        #time_full = stop_full - start_full

        #start_sm = time.time()
        li_gauss = prob.gaussian(self.const.METAL,li,measure_err)
        #li_gauss = norm.pdf(self.const.METAL,loc=li,scale=measure_err)
        mask = li_gauss > prob.FIVE_SIGMAS
        li_gauss = li_gauss[mask]
        li_gauss = li_gauss.reshape(len(li_gauss),1)
        #print "new shape",li_gauss.shape
        METAL = self.const.METAL[mask]
        METAL = METAL.reshape(len(METAL),1)

        #start = time.time()
        rv = lognorm(s=sig,scale=mu)
        astro_gauss = rv.pdf(METAL)
        product = li_gauss*astro_gauss
        #end = time.time()
        #print("Product Time",end - start)
        #print astro_gauss.shape
        
        #start = time.time()
        integral = prob.normalize(self.const.AGE,np.trapz(product,METAL,axis=0))
        #end = time.time()
        #print("Integral Time",end - start)
        #stop_sm = time.time()
        #time_small = stop_sm - start_sm

        return integral

    #calculates and returns a 2D array of sigma b-v and age
    #omit_cluster specifies a cluster index to remove from the fits to make the grids without a cluster
    def make_grids(self,bv_li,fits,upper_lim=None,medianSavefile=None,sigmaSavefile=None,\
            setAsDefaults=False, omit_cluster=None):
        if (medianSavefile and medianSavefile[-4:] == '.npy'):
            medianSavefile = medianSavefile[:-4]
        if (sigmaSavefile and sigmaSavefile[-4:] == '.npy'):
            sigmaSavefile = sigmaSavefile[:-4]

        primordial_li_fit = None
        bldb_fit = None
        li_scatter_fit = None #fit as function of LiEW
        ca_scatter = None #constant across all bv so done ahead of time
        fit_gaussian = None
        if (self.metal == 'lithium'):
            #primordial_li_fit = my_fits.primordial_li()
            primordial_li_fit = my_fits.MIST_primordial_li()#ngc2264_fit),fromFile=False,saveToFile=True)
            bldb_fit = my_fits.bldb_fit(fits)
            li_scatter_fit = my_fits.fit_two_scatters(bv_li,fits,upper_lim=upper_lim,\
                    omit_cluster=omit_cluster)
        elif (self.metal == 'calcium'):
            ca_scatter = my_fits.total_scatter(bv_li,fits,omit_cluster)
            clust_scatter = [fits[i][1](0.65) for i in range(len(fits))]
            fit_gaussian = my_fits.fit_gaussian(np.log10(self.const.CLUSTER_AGES),clust_scatter)
        
        median_rhk, sigma = [],[] #holds the grids
        for bv in self.const.BV_S:
            rhk,scatter,CLUSTER_AGES,_ = my_fits.get_valid_metal(bv,fits,self.const,\
                    primordial_li_fit,bldb_fit,omit_cluster)
            
            mu,sig,_ = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,self.metal,\
                    li_scatter_fit,ca_scatter,fit_gaussian)
            median_rhk.append(mu(self.const.AGE))
            sigma.append(sig(self.const.AGE))
            if self.metal=='calcium':        ########## SHORTCUT #############
                break
            
        median_rhk,sigma = np.array(median_rhk),np.array(sigma)

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


if  __name__ == "__main__":
    err = "Usage:  python2.7 baffles.py -bv <B-V> -rhk <Log10(R\'HK)> -li <Log10(EW/mA)> [-ul -maxAge <13000> -bverr <> -noplot -s -showplots -filename <> -help]"
    
    help_msg = "\n\
    -bv corrected (B-V)o of the star \n\
    -rhk <> log base 10 of the R'HK value \n\
    -li <> log base 10 of the EW measure in milli-angstroms (=0.1 pm) \n\
    \noptional flags:\n\
    -ul indicates that the log(EW/mA) reading is an upper-limit reading. \n\
    -maxAge allows user to input max posible age of star (Myr) if upper-limit flag is used default. is 13000 \n\
    -bverr <float uncertainity> provides the uncertainty on B-V with default .002 \n\
    -noplot suppresses all plotting and just generates .csv files (-s extension is implied. Used with -showplots it prevents saving).\n\
    -s saves posteriors as .csv \n\
    -showplots to show posterior plots to user before saving. \n\
    -filename <name of file w/o extension> name of files to be saved: name.pdf is graphs, name_calcium.csv/name_lithium.csv/name_product.csv are posterior csv files for calcium/lithium/product respectively.  csv is stored as age, probability in two 13000 element columns.\n\
    -help prints this message \n"
    
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
    upperLim = False
    maxAge = const.GALAXY_AGE
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
        if ('-showplots' in argv or '-showplot' in argv):
            showPlots = True
        if ('-noplot' in argv or '-noplots' in argv):
            noPlots = True
            save = True
        if ('-filename' in argv):
            fileName = argv[argv.index('-filename') + 1]
        if ('-bverr' in argv):
            bv_err = float(argv[argv.index('-bverr') + 1])
        if ('-ul' in argv or '-UL' in argv):
            upperLim = True
        if ('-maxAge' in argv):
            maxAge = float(argv[argv.index('-maxAge') + 1])
    except IndexError:
        print err
    except ValueError:
        print err
    
    baffles_age(bv,rhk,li,bv_err,upperLim,maxAge,fileName,showPlots=showPlots,noPlots=noPlots, savePostAsText=save)
