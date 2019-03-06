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

def main():
    #METAL = "lithium"
    METAL = "calcium"
    bv_m,upper_lim,fits = None,None,None
    if (METAL=='calcium'):
        global const
        import ca_constants as const
        bv_m,fits = readData.read_calcium(fromFile=True)#, saveToFile=False,fit_degree=0)
    else :
        bv_m, upper_lim, fits = readData.read_lithium(fromFile=True)


    #----------------------------------------------------------------------------------------------------
    #Plot Fits
    #my_plot.plot_fits(bv_m,fits,METAL,showPlots=True,upper_lim=upper_lim)
    
    #Metal vs B-V
    #pp=PdfPages(METAL + '_metal_vs_bv.pdf')
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True)
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True)
    #pp.close()

    #Metal vs Age
    #pp=PdfPages(METAL + '_metal_v_age.pdf')
    #bvRange = range(.24,2.3,1) if METAL = 'lithium' else range(.45,.9)
    #for bv in bvRange:
        #my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=True,shadeScatter=True,errorbars=True)
    #pp.close()

    #Scatter vs B-V
    #pp=PdfPages(METAL + '_scatter_vs_bv_with_offset.pdf')
    #my_plot.scatter_vs_bv(fits,METAL,pp,showPlots=True)
    #pp.close()

    #Scatter vs Age
    #pp=PdfPages(METAL + '_scatter_vs_age.pdf')
    #for bv in [.65]:
        #my_plot.scatter_vs_age(fits,METAL,bv,pp,showPlots=True)
    #pp.close()

    #Fit histogram
    #pp=PdfPages(METAL + '_hist_range.pdf')
    #for i in np.arange(1,3,.5):
    #    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=True,upper_limits=upper_lim,li_range=[i,i+.5],title="Log(LiEW) between %.1f and %.1f" % (i,i+.5))
    #pp.close()
    
    



    #Omitting each cluster
    pp = PdfPages(METAL + '_omit_clusters.pdf')
    for i in range(1,len(bv_m) - 1):
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        baf.posterior_product(bv_m[i][0],bv_m[i][1],pp,showPlot=True,showStars=False,givenAge=const.CLUSTER_AGES[i],title= const.CLUSTER_NAMES[i] + ' Posterior Product')
    pp.close()


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
