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
import ca_constants as const
import plotting as my_plot
import readData

def main():
    bv_rhk,fits = readData.read_calcium(fromFile=False)
    bv_li, upper_lim, li_fits = readData.read_lithium()
    bv_li_lin = []
    for c in bv_li:
        bv_li_lin.append([c[0],np.power(10,c[1])])
    li_fits = readData.get_li_fits(bv_li_lin,upper_lim)
    name = 'HD 984'
    bv = 7.82 - 7.32
    bv_err = np.sqrt(.02*.02 + .01*.01)
    rhk = -4.42
    li = 120 #np.log10(120)
    
    pp = PdfPages('hd984_rhk_li.pdf')
    plt.text(bv - .03,rhk + .09,name,size=13,color='r')
    plt.plot(bv,rhk,marker='*',markersize=18,color='r',linestyle='None',zorder=10)
    #my_plot.metal_vs_bv(bv_rhk,fits,'calcium',pp,True)
    my_plot.plot_mamajek(bv_rhk,fits)
    pp.savefig()
    plt.show()

    plt.text(bv - .11, li + 12,name,color='r',size=13)
    plt.plot(bv,li,marker='*',markersize=18,color='r',linestyle='None',zorder=10)
    my_plot.metal_vs_bv(bv_li_lin,li_fits,'lithium',pp,True,upper_lim)
    
    pp.close()
    

    """
    #creating posteriors
    #Equivalent to baffles.baffles_age(bv,rhk,li,name,pp,showPlots=True)
    pp = PdfPages('graphs/HD984_baffles.pdf')   
    baf = baffles.age_estimator('calcium')
    p = baf.get_posterior(bv,rhk,pp,True,bv_uncertainty=bv_err,mamajekAge=baffles.getMamaAge(rhk))
    #np.savetxt("hd984_calcium.csv", zip(const.AGE,p.array), delimiter=",")

    baf2 = baffles.age_estimator('lithium')
    p2 = baf2.get_posterior(bv,li,pp,True,bv_uncertainty=bv_err)
    #np.savetxt("hd984_lithium.csv", zip(const.AGE,p2.array), delimiter=",")

    title = name + ' Calcium/Lithium Posterior Product'
    y = p.array * p2.array
    prob.normalize(const.AGE,y)
    my_plot.posterior(const.AGE,y,prob.stats(const.AGE,y),title,pp,True)    
    #np.savetxt("hd984_posterior_product.csv", zip(const.AGE,y), delimiter=",")
    pp.close()
    """

if  __name__ == "__main__":
    main()
