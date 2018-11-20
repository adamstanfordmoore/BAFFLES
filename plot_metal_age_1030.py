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
    bv_li, upper_lim, li_fits = readData.read_lithium()
    #bv_rhk,rhk_fits = readData.read_calcium()

    
    pp = PdfPages('li_age_constrained_polynomial.pdf')
    for bv in np.arange(1,2,.1):
        my_plot.metal_vs_age(li_fits,'lithium',bv,pp,showPlots=True,title='Constrained Polynomial at B-V= %s' % bv)
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
