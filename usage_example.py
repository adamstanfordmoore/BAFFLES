"""
Adam Stanford-Moore
8/28/18
This file provides an example program for how to use the baffles module
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
    # Example 1: age of the sun
    bmv = 0.65 
    rhk = -4.906
        


    baf_ca = baffles.age_estimator('calcium')
    
    
    
    
    
    baf_li = baffles.age_estimator('lithium')

    






if  __name__ == "__main__":
    main()
