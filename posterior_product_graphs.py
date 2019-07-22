"""
Adam Stanford-Moore
8/28/18
This file provides an example program for how to use the baffles module
"""
import warnings
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import copy
import sys
import fitting as my_fits
import probability as prob
import baffles
import li_constants as const
import plotting as my_plot
import readData
from readData import in_bounds

"""
def in_bounds(bv,l,const,log=False):
    if (bv > const.BV_RANGE[0] and bv < const.BV_RANGE[1]):
        if (not log and (l > np.power(10,const.METAL_RANGE[0]) and l < np.power(10,const.METAL_RANGE[1]))):
            return True
        elif (log and (l > const.METAL_RANGE[0] and l < const.METAL_RANGE[1])):
            return True
        return False
"""

def main():
    bv_li, upper_lim, li_fits = readData.read_lithium()
    #bv_rhk,rhk_fits = readData.read_calcium()#fromFile=False, saveToFile=False)
    
    teff_to_bv = my_fits.magic_table_convert(1,6)
    def BV_ERR(t,err):
        return np.abs(teff_to_bv(t + err) - teff_to_bv(t-err))/2
    
    baf = baffles.age_estimator('lithium')#,'grids/median_li_071119','grids/sigma_li_071119')
    pp = PdfPages('max_bv_tuchor.pdf')
    for MAX_BV in [2]:# [1.7,1.4,1.3,1.1]:
        abdor_c = []
        abdor_l = []
        abdor_bverr = []
        abdor_lerr = []
        t = ascii.read('data/ab_dor_updated_err.csv',delimiter=',')
        #t = ascii.read('data/ab_dor.csv',delimiter=',')
        for line in t[1:]:
            bv = teff_to_bv(float(line[5]))
            ew = float(line[3])
            if in_bounds(bv,ew,const) and bv < MAX_BV:
                abdor_c.append(bv)
                abdor_l.append(ew)
                abdor_lerr.append(float(line[4]))
                abdor_bverr.append(BV_ERR(float(line[5]),float(line[6])))
                #lim_bp.append(False)
        abdor_c,abdor_l = np.array(abdor_c),np.log10(np.array(abdor_l))
        
        tuchor_c = []
        tuchor_l = []
        tuchor_bverr,tuchor_lerr = [],[]
        t = ascii.read('data/tuchor_updated_err.csv',delimiter=',')
        #t = ascii.read('data/tuchor.csv',delimiter=',')
        for line in t[1:]:
            bv = teff_to_bv(float(line[5]))
            ew = float(line[3])
            if in_bounds(bv,ew,const) and bv < MAX_BV:
                tuchor_c.append(bv)
                tuchor_l.append(ew)
                tuchor_lerr.append(float(line[4]))
                tuchor_bverr.append(BV_ERR(float(line[5]),float(line[6])))
                #lim_bp.append(False)
        tuchor_c,tuchor_l = np.array(tuchor_c),np.log10(np.array(tuchor_l))

        bp_c = []
        bp_l = []
        lim_bp = []
        bp_err = []
        bp_bv_err = []
        t = ascii.read('data/beta_pic_updated_err.txt', delimiter='\t')
        #t = ascii.read('data/beta_pic_noM.txt', delimiter='\t')
        for line in t:
            bv = teff_to_bv(float(line[7])) 
            #print "bv-diff:",np.abs(bv - float(line[4]))
            ew = float(line[5])
            if in_bounds(bv,ew,const) and ew != 20:
                bp_c.append(bv)
                bp_l.append(ew)
                #if ew== 20: 
                #    bp_l[-1] = 20
                #    lim_bp.append(True)
                #else: lim_bp.append(False)
                lim_bp.append(False)
                bp_err.append(float(line[6]))
                bp_bv_err.append(BV_ERR(float(line[7]),float(line[8])))
                #print "BV err",bp_bv_err[-1]
        bp_c,bp_l = np.array(bp_c),np.log10(np.array(bp_l))
        #bv_li.append([bp_c,bp_l])
        #upper_lim.append(lim_bp)
       
        #plt.scatter(bp_c,bp_l)
        #plt.plot(bp_c,my_fits.poly_fit(bp_c,bp_l)(bp_c))
        #plt.show()
        
        
        pp = PdfPages('moving_group_age.pdf')
        baf.posterior_product(abdor_c,abdor_l,abdor_bverr,abdor_lerr,pdfPage=pp,showPlot=True,showStars=True,title='AB DOR Posterior Product: MAX_BV= %.2f' % MAX_BV,givenAge=149,givenErr=[-19,51])
        baf.posterior_product(tuchor_c,tuchor_l,tuchor_bverr,tuchor_lerr,pdfPage=pp,showPlot=True,showStars=True,title='TUC/HOR Posterior Product: MAX_BV= %.2f' % MAX_BV,givenAge=45,givenErr=4)
        #pp.close()


        #baf = baffles.age_estimator('li',default_grids=False)#'grids/median_li_071119','grids/sigma_li_071119',default_grids=False)
        #baf.make_grids(bv_li,li_fits,omit_cluster=1)
        
        #mask = (bp_c < 1) & (bp_c > 0.6)
        #mask = (abdor_c < 1.2) & (abdor_c > 0.8)
        #for b,l in zip(bp_c[mask],bp_l[mask]):
        #for b,l in zip(abdor_c[mask],abdor_l[mask]):
        #    print b,l
            #my_plot.metal_vs_age(li_fits,'lithium',.72,pp,showPlots=True,shadeScatter=True,errorbars=True,title='B-V= %s' % .72, bv_m=bv_li,upper_lim=upper_lim,metal_val=l,omit_cluster=1)
            #my_plot.metal_vs_age(li_fits,'lithium',.704,pp,showPlots=True,shadeScatter=True,errorbars=True,title='B-V= %s' % .704, bv_m=bv_li,upper_lim=upper_lim,metal_val=l,omit_cluster=1)
            #my_plot.metal_vs_age(li_fits,'lithium',b,pp,showPlots=True,shadeScatter=True,errorbars=True,title='B-V= %s' % b, bv_m=bv_li,upper_lim=upper_lim,metal_val=l,omit_cluster=1)
        #    p = baf.get_posterior(b,l,pdfPage=pp,showPlot=True,logPlot=False,upperLim = False)
            




        #pp = PdfPages('beta_pic_posterior_product_v4.pdf')
        baf = baffles.age_estimator('li',default_grids=False)
        baf.make_grids(bv_li,li_fits,upper_lim,omit_cluster=1)
        #baf.resample_posterior_product(bp_c,bp_l,bp_bv_err,bp_err,lim_bp,pdfPage=pp,showPlot=True,showStars=True,givenAge=24,givenErr=3,title=r'$\beta$ Pic Posterior Product')
        baf.posterior_product(bp_c,bp_l,bp_bv_err,bp_err,lim_bp,pdfPage=pp,showPlot=True,showStars=True,givenAge=24,givenErr=3,title=r'$\beta$ Pic Posterior Product')
    pp.close()
        

if  __name__ == "__main__":
    main()
