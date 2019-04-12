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
import ca_constants as const
import plotting as my_plot
import readData

def main():
    bv_rhk,fits = readData.read_calcium()#fromFile=False)
    #bv_li, upper_lim, li_fits = readData.read_lithium()
    
    pp = PdfPages('baffles_vs_mamajek.pdf') 
    for i in range(len(bv_rhk)):
        my_plot.baffles_vs_mamajek(bv_rhk,fits,i,pp,showPlots=True,title=const.CLUSTER_NAMES[i])
    pp.close()

    """
    names = ['HD 984','HD 206893','HR 2562']
    bv = [0.5,.44,.45]
    rhk = [-4.42,-4.466,-4.551]
    
    pp = PdfPages('resultant_posteriors.pdf')
    baf = baffles.age_estimator('calcium')
    for i in range(len(names)):
        p = baf.get_posterior(bv[i],rhk[i],pp,True,mamajekAge=baffles.getMamaAge(rhk[i]),title=names[i]+ ' Calcium Posterior')
        plt.plot(bv[i],rhk[i],marker='*',markersize=18,color='r',linestyle='None',zorder=10)
        my_plot.metal_vs_bv(bv_rhk,fits,'calcium',pp,True)
    pp.close() 
    """
    

if  __name__ == "__main__":
    main()
