"""
Adam Stanford-Moore
2/11/19
"""
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import fitting as my_fits
import probability as prob
import baffles
import ca_constants as const
import plotting as my_plot
import readData
import utils
METAL = "calcium"
upper_lim = None
bv_m,fits = readData.read_calcium()#fromFile=False,saveToFile=False,fit_degree=0)
print([len(x[0]) for x in bv_m])

def main():
    plot_fits()
    metal_vs_bv()
    metal_vs_age()
    scatter_vs_age()
    fit_hist()
    #baffles_vs_mamajek()
    #combined_validation()




def printName(n):
    print("scp sasm@gpicruncher.stanford.edu:~/BAFFLES/" + n + "  ~/Desktop/")
    #print("sips -s format png -s formatOptions best ~/Desktop/" + n[6:] + " --out ~/Desktop/"

def plot_fits(cluster_indices=None):
    #Plot Fits
    name = 'plots/' + METAL + '_fits.pdf'
    pp=PdfPages(name)
    my_plot.plot_fits(bv_m,fits,METAL,pdfPage=pp,showPlots=False,upper_lim=upper_lim,specific_clusters=cluster_indices)
    pp.close()
    printName(name)

def metal_vs_bv():
    #Metal vs B-V
    name = 'plots/' +METAL + '_metal_vs_bv.pdf'
    pp=PdfPages(name)
    ##clusters = [0,4,6,7]
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True)#,specific_clusters = clusters)
    #my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True,shadeScatter=False)
    
    bv_m_linear,fits_linear = readData.read_calcium(fromFile=False,saveToFile=False,fit_degree=1)
    my_plot.metal_vs_bv(bv_m_linear,fits_linear,METAL,pp,showPlots=True)#,specific_clusters = clusters)
    #my_plot.metal_vs_bv(bv_m_linear,fits_linear,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True,shadeScatter=False)
    printName(name)
    pp.close()

def metal_vs_age():
    #Metal vs Age
    name = 'plots/' + METAL + '_metal_v_age.pdf'
    pp=PdfPages(name)
    bvRange = const.BV_S 
    for bv in bvRange:
        print("bv: ",bv)
        my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=True,shadeScatter=False,\
                errorbars=True,title=' ', bv_m=bv_m,upper_lim=upper_lim,\
                logAge=True,plotStars=False,mamajek_poly=True)
    printName(name)
    pp.close()

def scatter_vs_bv():
    #Scatter vs B-V
    pp=PdfPages('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    my_plot.scatter_vs_bv(fits,METAL,pp,showPlots=True)
    names.append('plots/' +METAL + '_scatter_vs_bv_with_offset.pdf')
    pp.close()

def scatter_vs_age():
    #Scatter vs Age
    name = 'plots/' +METAL + '_scatter_vs_age.pdf'
    pp=PdfPages(name)
    my_plot.scatter_vs_age(fits,METAL,.65,pp,showPlots=True,bv_m=bv_m,upper_lim=upper_lim)
    pp.close()
    printName(name)

def fit_hist():
    #Fit histogram
    name = 'plots/' +METAL + '_hist_fit.pdf'
    pp=PdfPages(name)
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=True)
    printName(name)
    pp.close()

def baffles_vs_mamajek():
    # BAFFLES vs Mamajek
    name = 'plots/baffles_vs_mamajek.pdf'
    pp = PdfPages(name)
    for i in range(len(bv_m)):
        my_plot.baffles_vs_mamajek(bv_m,fits,i,pdfPage=pp,showPlots=False,title=None,mamaProduct=False)
    pp.close()
    printName(name)

def omitting(validation=False):
    #Omitting each cluster
    #name = 'plots/' +METAL + '_self_validation.pdf'
    name = 'plots/' +METAL + '_omit_clusters.pdf' if not validation else \
            'plots/' +METAL + '_self_validation.pdf'
    pp = PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=validation)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        #mamaProductAge = utils.getMamaProductAge(bv_m[i][1])
        #print(mamaProductAge
        if not validation: baf.make_grids(bv_m,fits,omit_cluster=i)
        p = baf.posterior_product(bv_m[i][0],bv_m[i][1],pdfPage=pp,showPlot=True,\
                showStars=True,givenAge=const.CLUSTER_AGES[i],\
                title= const.CLUSTER_NAMES[i])
        #print("True: %.3g BAFFLES: %.3g Mama: %.3g" % (const.CLUSTER_AGES[i],p.stats[2],mamaProductAge)
    printName(name)
    pp.close()


def combined_validation():
    #Omitting each cluster
    name = 'plots/' +METAL + '_combined_validation.pdf'
    pp=PdfPages(name)
    baf_default = baffles.age_estimator(METAL)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p_val = baf.posterior_product(bv_m[i][0],bv_m[i][1])

        plt.plot(const.AGE,p_val.array,linewidth=2,linestyle='--',label='Posterior with removal')
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],\
            pdfPage=pp,showPlot=False,\
            showStars=True,givenAge=const.CLUSTER_AGES[i],\
            title= const.CLUSTER_NAMES[i])
        print(np.sum(p_val.array - p.array))
    printName(name)
    pp.close()



def posteriors():
    #Making posteriors
    name = 'plots/' + METAL + '_posteriors_tail.pdf'
    pp = PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=False)
    baf.make_grids(bv_m,fits,upper_lim)#,omit_cluster=0)
    #for bv in [0.65]:#[1.5,1.75,1.9]:
    #    for li in np.linspace(-3.8,-5,5):#[1,2,3]:#[3.1,1]:
    for bv,li in zip(bv_m[0][0],bv_m[1][1]):
        #my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=True,shadeScatter=False,errorbars=True,title='B-V= %s' % bv, bv_m=bv_m,upper_lim=upper_lim,metal_val=li)
    #        #baffles.baffles_age(bv,li=li,showPlots=True,pdfPage=pp)
        p = baf.get_posterior(bv,li,pdfPage=pp,showPlot=True,logPlot=False,upperLim = False,mamajekAge=True)
    printName(name)
    pp.close()
    





   
if  __name__ == "__main__":
    main()
