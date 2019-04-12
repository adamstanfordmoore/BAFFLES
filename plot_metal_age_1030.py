"""
Adam Stanford-Moore
8/28/18
This file provides an example program for how to use the baffles module
"""
import warnings
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
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
    bv_li, upper_lim, li_fits = readData.read_lithium(fromFile=False)
    bv_rhk,rhk_fits = readData.read_calcium(fromFile=False)

    
    #pp = PdfPages('li_scatter_v_age.pdf')
    #for bv in np.arange(.5,2,.5):
    #    my_plot.scatter_vs_age(li_fits,'lithium',bv,pp,showPlots=False,title='B-V= %s' % bv)
    #pp.close() 
    """
    pp = PdfPages('li_age_fits_all.pdf')
    for bv in np.arange(1,2,.05):
        my_plot.metal_vs_age(li_fits,'lithium',bv,pp,showPlots=True,title='B-V= %s' % bv,errorbars=True)
    pp.close() 
    """ 
    """
    pp = PdfPages('postBLDB_posteriors')
    for bv in [1.5,1.75,2,2.25]:
        for li in [2,1.2]:
            baffles.baffles_age(bv,li=li,pdfPage=pp,showPlots=True)
    pp.close()
    """ 

    
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,rhk_fits,'grids/median_rhk_030719','grids/sigma_rhk_030719',True)
    
    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(bv_li,li_fits,'grids/median_li_030719','grids/sigma_li_030719',True,upper_lim=upper_lim)
    

if  __name__ == "__main__":
    main()
