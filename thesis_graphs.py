"""
Adam Stanford-Moore
2/11/19
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
    #bv_li, upper_lim, li_fits = readData.read_lithium(fromFile=True)
    bv_rhk,rhk_fits = readData.read_calcium(fromFile=False, saveToFile=False,fit_degree=0)
    
    #Scatter vs B-V
    #pp=PdfPages('scatter_vs_bv_with_offset.pdf')
    #my_plot.scatter_vs_bv(li_fits,'li',pp,showPlots=True)
    #pp.close()


    #Fit histogram
    pp=PdfPages('ca_constant_hist.pdf')
    my_plot.fit_histogram(bv_rhk,rhk_fits,'ca',pp,showPlots=True)
    pp.close()


    #pp=PdfPages('poster_metal_v_age.pdf')
    #my_plot.metal_vs_age(rhk_fits,'ca',.65,pp,showPlots=True,shadeScatter=True,errorbars=True)
    #pp.close()


    #pp = PdfPages('mamajek_calcium.pdf')
    #my_plot.metal_vs_bv(bv_rhk,rhk_fits,'calcium',pp,showPlots=True)
    #pp.close() 
    """
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

    pp = PdfPages('beta_pic_posterior_product_v2.pdf')
    baf = baffles.age_estimator('li')
    baf.make_grids(li_fits)
    baf.posterior_product(bp_c,bp_l,pp,showPlot=True,showStars=True,givenAge=24,title=' ')
    pp.close()
    """
    """ 
    #Omitting each cluster
    pp = PdfPages('clusters_5clusters.pdf')
    for i in range(1,len(bv_li) - 1):
        print "Omitting " + const.CLUSTER_NAMES[i]
        baf = baffles.age_estimator('li',default_grids=False)
        baf.make_grids(li_fits,omit_cluster=-1)
        baf.posterior_product(bv_li[i][0],bv_li[i][1],pp,showPlot=True,showStars=False,givenAge=const.CLUSTER_AGES[i],title= const.CLUSTER_NAMES[i] + ' Posterior Product')
    pp.close()
    """


    """
    pp = PdfPages('postBLDB_posteriors')
    for bv in [1.5,1.75,2,2.25]:
        for li in [2,1.2]:
            baffles.baffles_age(bv,li=li,pdfPage=pp,showPlots=True)
    pp.close()
    """ 

    """
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(rhk_fits,'grids/median_rhk_103018','grids/sigma_rhk_103018',True)

    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(li_fits,'grids/median_li_103018','grids/sigma_li_103018',True)
    """

if  __name__ == "__main__":
    main()
