"""
Adam Stanford-Moore
8/15/19
This file provides examples for how to use the baffles module
"""
import numpy as np
import matplotlib.pyplot as plt
import fitting as my_fits
import probability as prob
import baffles
import plotting as my_plot
import readData
from matplotlib.backends.backend_pdf import PdfPages

def main():
    # Example 1: Single line to find the age of the sun
    baffles.baffles_age(bv=0.65,rhk=-4.906)

    # Example 2: Longer line to find age of HD 206893 and save
    baffles.baffles_age(bv=0.45,rhk=-4.55,li=21,bv_err=.02,li_err = 5,
                        upperLim=False, maxAge=None,fileName='HD206893_example',
                        pdfPage=None, showPlots=True,noPlots=False,
                        savePostAsText=False)


    # Example 3: Create posterior and use posterior after
    # first load baffles class for calcium which contains functions to find posterior
    baffles_ca = baffles.age_estimator('calcium')

    # posterior is a struct with array, stats (array of ages at CDF .02,.16,.5,.84,.97)
    posterior = baffles_ca.get_posterior(bv=0.65,metallicity=-4.906)

    # can manipulate posterior.array as you would like
    plt.plot(ca_const.AGE,posterior.array)
    plt.show()
    # or use built-in plotting function
    my_plot.posterior(ca_const.AGE,posterior.array,posterior.stats,showPlot=True)




    # Exmple 3: Lets compute a posterior product and save it as a pdf
    tuchor_bmv, tuchor_li,tuchor_li_err = readData.tuchor()
    baffles_li = baffles.age_estimator('lithium')
    pp = PdfPages("tuchor_example_product.pdf")
    product = baffles_li.posterior_product(tuchor_bmv,tuchor_li,bv_errs=None,
                                 tuchor_li_err,pdfPage=pp,showPlot=True,
                                 showStars=True,title='Tuc/Hor',
                                 givenAge=45,givenErr=4)
    # product.array contains posterior product as before


    
    
    
    
    


    






if  __name__ == "__main__":
    main()
