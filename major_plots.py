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
import fitting as my_fits

def main():
    #METAL = "lithium"
    METAL = "calcium"
    names,pp = [],None
    bv_m,upper_lim,fits = None,None,None
    if (METAL=='calcium'):
        global const
        import ca_constants as const
        bv_m,fits = readData.read_calcium()#, saveToFile=False,fit_degree=0)
    else :
        bv_m, upper_lim, fits = readData.read_lithium()#fromFile=False,saveToFile=True)

    print METAL
    for c in range(len(bv_m)):
        print const.CLUSTER_NAMES[c], len(bv_m[c][0])

    #----------------------------------------------------------------------------------------------------
    #Plot Fits
    #pp=PdfPages('plots/' + METAL + '_fits.pdf')
    #my_plot.plot_fits(bv_m,fits,METAL,showPlots=True,upper_lim=upper_lim)
    #names.append('plots/' +METAL + '_fits.pdf')
    #pp.close()

    #Metal vs B-V
    #pp=PdfPages('plots/' +METAL + '_metal_vs_bv.pdf')
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True)
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True)
    #names.append('plots/' +METAL + '_metal_vs_bv.pdf')
    #pp.close()

    #Metal vs Age
    #pp=PdfPages('plots/' +METAL + '_metal_v_age_041019.pdf')
    #bvRange = np.arange(.35,1.9,.1) if METAL == 'lithium' else np.arange(.45,.9,.05) #np.arange(.35,1.9,.1)
    #for bv in bvRange:
    #    my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=False,shadeScatter=True,errorbars=True,title='B-V= %s' % bv, from_fit_two_scatters=True,bv_m=bv_m,upper_lim=upper_lim)
    #names.append('plots/' +METAL + '_metal_v_age_041019.pdf')
    #pp.close()

    #Scatter vs B-V
    #pp=PdfPages('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    #my_plot.scatter_vs_bv(fits,METAL,pp,showPlots=True)
    #names.append('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    #pp.close()

    #Scatter vs Age
    #pp=PdfPages('plots/' +METAL + '_scatter_vs_age.pdf')
    #for bv in [.65]:
    #    my_plot.scatter_vs_age(fits,METAL,bv,pp,showPlots=True)
    #names.append('plots/' +METAL + '_scatter_vs_age.pdf')
    #pp.close()

    #Fit histogram
    #pp=PdfPages('plots/' +METAL + '_hist_range.pdf')
    #for i in np.arange(1,3,.5):
    #    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=True,upper_limits=upper_lim,li_range=[i,i+.5],title="Log(LiEW) between %.1f and %.1f" % (i,i+.5))
    #names.append('plots/' +METAL + '_hist_range.pdf')
    #pp.close()
    
    #BLDB fit
    #my_fits.bldb_fit(fits,plot=True)
    #names.append('bldb_vs_age.pdf')
    

    """
    #Omitting each cluster
    pp = PdfPages('plots/' +METAL + '_omit_clusters.pdf')
    for i in range(1,len(bv_m) - 1):
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        baf.posterior_product(bv_m[i][0],bv_m[i][1],pp,showPlot=False,showStars=False,givenAge=const.CLUSTER_AGES[i],title= const.CLUSTER_NAMES[i] + ' Posterior Product')
    #names.append('plots/' +METAL + '_omit_clusters.pdf')
    pp.close()
    """
    """
    pp = PdfPages(''plots/' +postBLDB_posteriors')
    for bv in [1.5,1.75,2,2.25]:
        for li in [2,1.2]:
            baffles.baffles_age(bv,li=li,pdfPage=pp,showPlots=True)
    #names.append(''plots/' +postBLDB_posteriors')
    pp.close()
    """ 

    """
    # Update default grids
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,rhk_fits,'grids/median_rhk_030719','grids/sigma_rhk_030719',True)

    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(bv_li,li_fits,'grids/median_li_030719','grids/sigma_li_030719',True,upper_lim=upper_lim)
    """

    for n in names:
        print "scp sasm@gpicruncher.stanford.edu:~/BAFFLES/" + n + "  ~/Desktop/"
        print "sips -s format png -s formatOptions best ~/Desktop/" + n[6:] + " --out ~/Desktop/"

if  __name__ == "__main__":
    main()
