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
import utils

def main():
    #METAL = "lithium"
    METAL = "calcium"
    names,pp = [],None
    bv_m,upper_lim,fits = None,None,None
    if (METAL=='calcium'):
        global const
        import ca_constants as const
        bv_m,fits = readData.read_calcium(fromFile=False,saveToFile=False)#, saveToFile=False,fit_degree=0)
    else :
        bv_m, upper_lim, fits = readData.read_lithium()#fromFile=False,saveToFile=True)
        #bv_m, upper_lim, fits = bv_m[-8:], upper_lim[-8:], fits[-8:]

    #----------------------------------------------------------------------------------------------------
    
    #Plot Fits
    #name = 'plots/' + METAL + '_fits.pdf'
    #pp=PdfPages(name)
    #my_plot.plot_fits(bv_m,fits,METAL,pdfPage=pp,showPlots=True,upper_lim=upper_lim)#,specific_clusters=[9])
    #names.append(name)
    #pp.close()

    #Metal vs B-V
    #pp=PdfPages('plots/' +METAL + '_metal_vs_bv.pdf')
    #clusters = [0,2,4,6,8,9] if METAL == 'lithium' else [0,4,6,7]
    ##sets = [clusters[0:i] for i in range(1,len(clusters)+1)]
    ##for clus in sets:
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True)
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True,shadeScatter=False)
    #names.append('plots/' +METAL + '_new_prim_li.pdf')
    #pp.close()

    #Metal vs Age
    #name = 'plots/' + METAL + '_metal_v_age_0718.pdf'
    #pp=PdfPages(name)
    #bvRange = const.BV_S #np.linspace(0.3,1.9,40) #np.arange(1.3,2,.1)#np.arange(0.3,2,.1) if METAL == 'lithium' else np.arange(.45,.9,.05) #np.arange(.35,1.9,.1)
    #for bv in bvRange:
    #    print "bv: ",bv
    #    my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=True,shadeScatter=False,\
    #            errorbars=True,title='B-V= %s' % bv, bv_m=bv_m,upper_lim=upper_lim,\
    #            logAge=True,plotStars=True)
    #names.append(name)
    #pp.close()

    #Scatter vs B-V
    #pp=PdfPages('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    #my_plot.scatter_vs_bv(fits,METAL,pp,showPlots=True)
    #names.append('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    #pp.close()

    #Scatter vs Age
    #pp=PdfPages('plots/' +METAL + '_scatter_vs_age.pdf')
    #for bv in [1.4]:#[.65]:
    #    my_plot.scatter_vs_age(fits,METAL,bv,pp,showPlots=True,bv_m=bv_m,upper_lim=upper_lim)
    #names.append('plots/' +METAL + '_scatter_vs_age.pdf')
    #pp.close()

    #Fit histogram
    #pp=PdfPages('plots/' +METAL + '_hist_range.pdf')
    #for i in np.arange(1,3,.5):
    #my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=True,upper_limits=upper_lim)#, li_range=[np.log10(30),3.5])#,title="Log(LiEW) between %.1f and %.1f" % (i,i+.5))
    #names.append('plots/' +METAL + '_hist_range.pdf')
    #pp.close()
    
    #BLDB fit
    #my_fits.bldb_fit(fits,plot=True)
    #names.append('bldb_vs_age.pdf')
   
    # BAFFLES vs Mamajek
    #pp = PdfPages('plots/' +METAL + '_baffles_vs_mamajek.pdf')
    #for i in range(len(bv_m)):
    #    my_plot.baffles_vs_mamajek(bv_m,fits,i,pdfPage=pp,showPlots=False,title=None,mamaProduct=True)
    #pp.close()
    #names.append('plots/' +METAL + '_baffles_vs_mamajek.pdf')

    #Omitting each cluster
    #name = 'plots/' +METAL + '_self_validation.pdf'
    name = 'plots/' +METAL + '_omit_clusters_072219.pdf'
    pp = PdfPages(name)
    for i in range(len(bv_m)):
        print const.CLUSTER_NAMES[i]
    #    #mamaProductAge = utils.getMamaProductAge(bv_m[i][1])
    #    #print mamaProductAge
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p = baf.posterior_product(bv_m[i][0],bv_m[i][1],pdfPage=pp,showPlot=True,\
                showStars=True,givenAge=const.CLUSTER_AGES[i],\
                title= const.CLUSTER_NAMES[i] + ' Posterior Product')
    #    #print "True: %.3g BAFFLES: %.3g Mama: %.3g" % (const.CLUSTER_AGES[i],p.stats[2],mamaProductAge)
    names.append(name)
    pp.close()
    
    # Making posteriors
    #name = 'plots/' + METAL + '_posteriors_tail.pdf'
    #pp = PdfPages(name)
    #baf = baffles.age_estimator(METAL,default_grids=False)#,'grids/median_'+METAL[:2]+'_071119','grids/sigma_'+METAL[:2]+'_071119') #'grids/median_rhk_062719','grids/sigma_rhk_062719',default_grids=False)
    #baf.make_grids(bv_m,fits,upper_lim,'grids/median_'+METAL[:2]+'_072219','grids/sigma_'+METAL[:2]+'_072219',setAsDefaults=True)
    #for bv in [0.65]:#[1.5,1.75,1.9]:
    #    for li in np.linspace(-3.8,-5,5):#[1,2,3]:#[3.1,1]:
    #        my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=False,shadeScatter=False,errorbars=True,title='B-V= %s' % bv, bv_m=bv_m,upper_lim=upper_lim,metal_val=li)
    #        #baffles.baffles_age(bv,li=li,showPlots=True,pdfPage=pp)
    #        p = baf.get_posterior(bv,li,pdfPage=pp,showPlot=True,logPlot=False,upperLim = False,mamajekAge=True)
    #names.append(name)
    #pp.close()
    

    """
    # Update default grids
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,rhk_fits,'grids/median_rhk_030719','grids/sigma_rhk_030719',True)
    
    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(bv_m,fits,upper_lim,'grids/median_li_071819','grids/sigma_li_071819',setAsDefaults=True)
    """

    for n in names:
        print "scp sasm@gpicruncher.stanford.edu:~/BAFFLES/" + n + "  ~/Desktop/"
        #print "sips -s format png -s formatOptions best ~/Desktop/" + n[6:] + " --out ~/Desktop/"

if  __name__ == "__main__":
    main()
