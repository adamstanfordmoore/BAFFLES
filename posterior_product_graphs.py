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

def in_bounds(bv,l,const,log=False):
    if (bv > const.BV_RANGE[0] and bv < const.BV_RANGE[1]):
        if (not log and (l > np.power(10,const.METAL_RANGE[0]) and l < np.power(10,const.METAL_RANGE[1]))):
            return True
        elif (log and (l > const.METAL_RANGE[0] and l < const.METAL_RANGE[1])):
            return True
        return False

def main():
    bv_li, upper_lim, li_fits = readData.read_lithium()
    #bv_rhk,rhk_fits = readData.read_calcium()#fromFile=False, saveToFile=False)
    
    teff_to_bv = my_fits.magic_table_convert(1,6)
    abdor_c = []
    abdor_l = []
    t = ascii.read('data/ab_dor.csv',delimiter=',')
    for line in t[1:]:
        bv = teff_to_bv(float(line[4]))
        ew = float(line[3])
        if in_bounds(bv,ew,const):
            abdor_c.append(bv)
            abdor_l.append(ew)
            #lim_bp.append(False)
    abdor_c,abdor_l = np.array(abdor_c),np.log10(np.array(abdor_l))
    
    tuchor_c = []
    tuchor_l = []
    t = ascii.read('data/tuchor.csv',delimiter=',')
    for line in t[1:]:
        bv = teff_to_bv(float(line[4]))
        ew = float(line[3])
        if in_bounds(bv,ew,const):
            tuchor_c.append(bv)
            tuchor_l.append(ew)
            #lim_bp.append(False)
    tuchor_c,tuchor_l = np.array(tuchor_c),np.log10(np.array(tuchor_l))

    pp = PdfPages('moving_group_age.pdf')
    baf = baffles.age_estimator('li')
    baf.posterior_product(abdor_c,abdor_l,pp,showPlot=True,showStars=True,title='AB DOR Posterior Product',givenAge=149,givenErr=[-19,51])
    baf.posterior_product(tuchor_c,tuchor_l,pp,showPlot=True,showStars=True,title='TUC/HOR Posterior Product',givenAge=45,givenErr=4)
    #pp.close()


    
    bp_c = []
    bp_l = []
    lim_bp = []
    t = ascii.read('data/beta_pic_noM.txt', delimiter='\t')
    for line in t:
        if in_bounds(float(line[4]),float(line[5]),const):
            bp_c.append(float(line[4]))
            bp_l.append(float(line[5]))
            lim_bp.append(False)
    bp_c,bp_l = np.array(bp_c),np.log10(np.array(bp_l))
    #bv_li.append([bp_c,bp_l])
    #upper_lim.append(lim_bp)
   
    #plt.scatter(bp_c,bp_l)
    #plt.plot(bp_c,my_fits.poly_fit(bp_c,bp_l)(bp_c))
    #plt.show()

    #pp = PdfPages('beta_pic_posterior_product_v4.pdf')
    baf = baffles.age_estimator('li')
    baf.make_grids(bv_li,li_fits,omit_cluster=1)
    baf.posterior_product(bp_c,bp_l,pp,showPlot=True,showStars=True,givenAge=24,givenErr=3,title=r'$\beta$ Pic Posterior Product')
    pp.close()
    

if  __name__ == "__main__":
    main()
