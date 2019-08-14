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
    upper_lim = None
    METAL = 'calcium'
    #METAL = 'lithium'
    if METAL == 'calcium': bv_m,fits = readData.read_calcium()#fromFile=False)
    else: bv_m, upper_lim, fits = readData.read_lithium()
     
    """ 
    pp = PdfPages('baffles_vs_mamajek.pdf') 
    for i in range(len(bv_rhk)):
        my_plot.baffles_vs_mamajek(bv_rhk,fits,i,pp,showPlots=True,title=const.CLUSTER_NAMES[i])
    pp.close()
    """
    names = ['HD 984','HD 206893','HR 2562']# ['HIP 23418']#
    bv = [0.5,.44,.45 ]#[1.4]
    rhk = [-4.42,-4.466,-4.551] #[np.log10(20)]
    

    pp = PdfPages('hd984_on_mamajek.pdf')
    #pp = PdfPages('sun.pdf')
    #pp = PdfPages('plots/upper_lim_test.pdf')
    baf = baffles.age_estimator(METAL)
    for i in range(len(names)):
        p = baf.get_posterior(None,rhk[i],upperLim=False,maxAge=None,pdfPage=pp,showPlot=True,title=names[i] + " Calcium Posterior" ,logPlot=False,mamajekAge=True)
        plt.plot(bv[i],rhk[i],marker='s',markersize=26,color='w',linestyle='None',zorder=9)
        plt.plot(bv[i],rhk[i],marker='*',markersize=25,color='r',linestyle='None',zorder=10)
        my_plot.metal_vs_bv(bv_m,fits,METAL,pp,True,specific_clusters=[0,4,6,7])
    pp.close() 
    print "scpgpi hd984_on_mamajek.pdf"
    
    """ 
    pp =PdfPages("likelihood_comparison_hd984.pdf")
    baf = baffles.age_estimator('calcium')
    p = baf.get_posterior(0.5,-4.42,showPlot=False,pdfPage=pp,title='',mamajekAge=baffles.getMamaAge(-4.42))
    my_plot.posterior(const.AGE,p.array,p.stats,logPlot=True,pp=pp,showPlot=True,title='',mamajekAge=baffles.getMamaAge(-4.42))
    plt.show()
    """
    



    #pp.close()



if  __name__ == "__main__":
    main()
