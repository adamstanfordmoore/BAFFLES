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
import li_constants as const
import plotting as my_plot
import readData
import fitting as my_fits
import utils
import bisect

#from readData import schkolnik_betaPic
METAL = "lithium"
bv_m, upper_lim, fits = readData.read_lithium()#fromFile=False,saveToFile=False)
bv_ca, ca_fits = readData.read_calcium()#fromFile=False,saveToFile=False)

def printName(n):
    print("scp sasm@gpicruncher.stanford.edu:~/BAFFLES/" + n + "  ~/Desktop/")
    #print "sips -s format png -s formatOptions best ~/Desktop/" + n[6:] + " --out ~/Desktop/"

  
def main():
    #baf = baffles.age_estimator(METAL)
    #pp=PdfPages('plots/lithium_grid_mean')
    #for i in range(0,1000,200):
    #    plt.imshow(baf.grid_median[:,i:i+200])
    #    plt.xlabel("%d < Age < %d" % (const.AGE[i],const.AGE[i+199]))
    #    plt.ylabel("%.2f < B-V < %.2f" % (.35,1.9))
    #    plt.colorbar()
    #    pp.savefig()
    #    plt.show()
    #pp.close()

    #notable_stars()
    metal_vs_age()
    #plot_fits()
    #metal_vs_bv()
    #combined_validation()
    #omitting()
    #omitting(validation=True)
    #moving_group()


def moving_group():
    name = 'plots/moving_group_age.pdf'
    pp = PdfPages(name)
    baf = baffles.age_estimator('lithium')
    abdor = readData.abdor()
    fit = my_fits.poly_fit(abdor[0],abdor[1],2,scatter=True)
    plt.title("AB Dor")
    #my_plot.plot_fits([abdor],[fit],METAL,pdfPage=pp,showPlots=False)

    baf.posterior_product(abdor[0],abdor[1],None,abdor[2],pdfPage=pp,showPlot=False,showStars=True,title='AB Dor',givenAge=149,givenErr=[-19,51])

    tuc = readData.tuchor()
    fit = my_fits.poly_fit(tuc[0],tuc[1],2,scatter=True)
    plt.title("Tuc/Hor")
    #my_plot.plot_fits([tuc],[fit],METAL,pdfPage=pp,showPlots=False)
    baf.posterior_product(tuc[0],tuc[1],None,tuc[2],pdfPage=pp,showPlot=False,showStars=True,title='Tuc/Hor',givenAge=45,givenErr=4)

    #baf.make_grids(bv_m,fits,upper_lim,omit_cluster=1)
    #bp_b,bp_l,bp_ul,bp_err,names = readData.merged_betaPic()
    #print len(bp_b)
    #for name,b,l,ul in zip(names,bp_b,bp_l,bp_ul):
    #    print name,'\t',b,'\t',l,'\t',ul
    #fit = my_fits.poly_fit(bp_b,bp_l,2,bp_ul)
    #plt.title(r"$\beta$ Pic")
    #my_plot.plot_fits([[bp_b,bp_l]],[fit],METAL,pdfPage=pp,showPlots=False,upper_lim=[bp_ul])
    #baf.posterior_product(bp_b,bp_l,None,bp_err,bp_ul,pdfPage=pp,showPlot=False,showStars=True,givenAge=24,givenErr=3,title=r'$\beta$ Pic')

    pp.close()
    printName(name)



def plot_fits():
    #Plot Fits
    name = 'plots/' + METAL + '_fits.pdf'
    pp=PdfPages(name)
    my_plot.plot_fits(bv_m,fits,METAL,pdfPage=pp,showPlots=True,upper_lim=upper_lim,shadeScatter=False,specific_clusters=[8])
    printName(name)
    pp.close()

def metal_vs_bv():
    #Metal vs B-V
    name = 'plots/' +METAL + '_metal_vs_bv.pdf'
    pp=PdfPages(name)
    ###clusters = [0,4,6,7]
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,specific_clusters = None,legend=False)
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=True,upper_lim=upper_lim,primordial_li=True,fits_only=True,shadeScatter=False)
    printName(name)
    pp.close()

def metal_vs_age():
    #Metal vs Age
    name = 'plots/' + METAL + '_metal_v_age.pdf'
    pp=PdfPages(name)
    bvRange = const.BV_S #np.linspace(0.3,1.9,40) #np.arange(1.3,2,.1)#np.arange(0.3,2,.1) if METAL == 'lithium' else np.arange(.45,.9,.05) #np.arange(.35,1.9,.1)
    for bv in bvRange:
        print("bv: ",bv)
        my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=False,shadeScatter=False,\
                errorbars=False,title='B-V = %.3f' % bv, bv_m=bv_m,upper_lim=upper_lim,\
                logAge=True,plotStars=True)
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
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=True,upper_limits=upper_lim)
    printName(name)
    pp.close()
   
def bldb():
    #BLDB fit
    my_fits.bldb_fit(fits,plot=True)
    printName('bldb_vs_age.pdf')


def combined_validation():
    #Omitting each cluster
    name = 'plots/' +METAL + '_combined_validation.pdf'
    pp=PdfPages(name)
    baf_default = baffles.age_estimator(METAL)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p_val = baf.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=upper_lim[i])

        plt.plot(const.AGE,p_val.array,linewidth=2,linestyle='--',label='Posterior with removal')
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=upper_lim[i],\
                pdfPage=pp,showPlot=False,\
                showStars=True,givenAge=const.CLUSTER_AGES[i],\
                title= const.CLUSTER_NAMES[i])
        
    printName(name)
    pp.close()



def omitting(validation=False):
    #Omitting each cluster
    #name = 'plots/' +METAL + '_self_validation.pdf'
    name = 'plots/' +METAL + '_omit_clusters.pdf' if not validation else \
           'plots/' +METAL + '_self_validation.pdf'
    pp=PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=False)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        if not validation: baf.make_grids(bv_m,fits,omit_cluster=i)
        p = baf.posterior_product(bv_m[i][0],bv_m[i][1],pdfPage=pp,showPlot=False,\
                showStars=True,givenAge=const.CLUSTER_AGES[i],\
                title= const.CLUSTER_NAMES[i])
    printName(name)
    pp.close()


def notable_stars():
    li_const = utils.init_constants('lithium')
    name = 'plots/notable_stars.pdf'
    names = ["HR 2562","HD 206893","TW PsA"]
    bv = [.45,.44,1.1]
    bv_err = [np.sqrt(.014**2+.01**2),np.sqrt(.02**2+.01**2),.03]
    rhk = [-4.55,-4.466,None]
    mamaAge = [utils.getMamaAge(rhk[0]),utils.getMamaAge(rhk[1]),None]
    li = [21,28.5,33]
    li_err = [5,7,2]
    markers = ['*','*','*']
    markerSize = 25
    colors = ['gold','green','darkmagenta']
    age = [None,None,440]
    age_range = [[300,900],[200,2100],[400,480]]
    pp=PdfPages(name)

    my_plot.metal_vs_bv(bv_ca,ca_fits,'calcium',None,False,specific_clusters=[0,5,7,8],legend=False)
    plt.plot(bv,rhk,marker='s',markersize=markerSize,color='w',linestyle='None',zorder=9)
    for i in [0,1]:
        plt.plot(bv[i],rhk[i],marker=markers[i],markersize=markerSize,color=colors[i],linestyle='None',zorder=10,label=names[i])
    plt.legend()
    plt.xlim([.42,.9])
    pp.savefig()
    #plt.show()
    plt.close()

    my_plot.metal_vs_bv(bv_m,fits,'lithium',None,False,upper_lim=upper_lim,specific_clusters=[0,1,4,6,8,9])#,legend=False)
    plt.plot(bv,li,marker='s',markersize=markerSize,color='w',linestyle='None',zorder=9)
    for i in [0,1,2]:
        plt.plot(bv[i],np.log10(li[i]),marker=markers[i],markersize=markerSize,color=colors[i],linestyle='None',zorder=10,label=names[i])
    plt.legend()
    pp.savefig()
    #plt.show()
    plt.close()

    baf_ca = baffles.age_estimator('calcium')
    baf_li = baffles.age_estimator('lithium')
    for i in [0,1]:
        plt.plot([],[],'C0',linewidth=2,label='Final Age Posterior')
    
        p_li = baf_li.get_posterior(bv[i],li[i],bv_uncertainty=bv_err[i],measure_err=li_err[i],upperLim=False)
        p_ca = baf_ca.get_posterior(None,rhk[i])
        product = prob.normalize(const.AGE,p_ca.array*p_li.array)
        my_plot.posterior(const.AGE, product, prob.stats(const.AGE,product),names[i],logPlot=False)
        plt.plot(const.AGE, p_ca.array,color='C3',label="Calcium Posterior")
        plt.plot(const.AGE, p_li.array,color='C2',label="Lithium Posterior")
        plt.axvspan(age_range[i][0],age_range[i][1], alpha=0.2, color='r',label=r'Literature age: %d - %d Myr' % tuple(age_range[i]),zorder=0)
        plt.axvline(x=mamaAge[i],color='C5',linestyle='--',label='MH08 age: %d' % mamaAge[i])
        #plt.xlim([0,490])
        plt.legend()
        pp.savefig()
        #plt.show()
        plt.close()

    
    plt.axvspan(age_range[-1][0],age_range[-1][1], alpha=0.2, color='r',zorder = 0)
    plt.axvspan(360-140,360+140, alpha=0.2, color='b',zorder=0)
    p_li = baf_li.get_posterior(bv[2],li[2],bv_uncertainty=bv_err[2],measure_err=li_err[2],upperLim=False)
    my_plot.posterior(const.AGE, p_li.array, p_li.stats,names[2],None,False, logPlot=False)
    plt.axvline(x=age[-1],color='r',label=r'Literature age: %d $\pm$ 40 Myr' % age[-1])
    plt.axvline(x=360,color='b',label=r"M'12 Li age: 360 $\pm$ 140 Myr")
    plt.xlim([-30,510])
    plt.legend()
    pp.savefig()
    #plt.show()
    plt.close()
    
    
    plt.axvline(x=age[-1],color='r',label=r'Literature age: %d $\pm$ 40 Myr' % age[-1])
    plt.axvspan(age_range[-1][0],age_range[-1][1], alpha=0.2, color='r',zorder = 0)
    robs_f = robs_fomalhaut()
    plt.plot(const.AGE,robs_f(const.AGE),'k--',label='Nielsen 2019 Fomalhaut PDF') 
    plt.plot(const.AGE,p_li.array,'g',label='BAFFLES Li posterior') 
   
    y = prob.normalize(const.AGE,p_li.array*robs_f(const.AGE))
    plt.plot(const.AGE,y,color = 'C0',linewidth=2,label='Final Age') 
    stat = prob.stats(const.AGE,y)
    plt.vlines(x=stat[2],ymin= 0,ymax= y[bisect.bisect_left(const.AGE,stat[2])], \
                    label='Final median age: %.3g Myr' % stat[2] ,color = 'orange')
    plt.fill_between(const.AGE,y, where= (const.AGE >= stat[1]) & (const.AGE <= stat[3]),color='.4', \
                        label='68%% CI: %.2g - %.2g' % (stat[1],stat[-2]))
    plt.title(names[2])
    plt.xlim([0,1200])
    plt.legend()
    pp.savefig()
    #plt.show()
    plt.close()


    



    printName(name)
    pp.close()

def robs_fomalhaut():
    import astropy
    from scipy.signal import savgol_filter
    t = astropy.io.fits.open('data/Fom-age-pdf.fits')
    data = t[0].data
    x = np.linspace(np.min(data)-10,np.max(data)+10,1000)
    cdf = np.array([(data < n).sum() for n in x],dtype='float')/len(x)
    cdf /= cdf[-1]

    smoothed = savgol_filter(cdf, 51, 3)

    pdf = np.gradient(smoothed)
    pdf = savgol_filter(pdf,101,3)
    prob.normalize(x,pdf)
    
    #plt.plot(x,pdf,label='PDF')
    #plt.hist(data,bins=100,density=True)
    #plt.show()

    pdf[0:2],pdf[-2:] = [0,0],[0,0]
    
    return my_fits.piecewise(x,pdf)







def posteriors():
    # Making posteriors
    name = 'plots/' + METAL + '_posteriors_tail.pdf'
    pp = PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=False)#,'grids/median_'+METAL[:2]+'_071119','grids/sigma_'+METAL[:2]+'_071119') #'grids/median_rhk_062719','grids/sigma_rhk_062719',default_grids=False)
    baf.make_grids(bv_m,fits,upper_lim,'grids/median_'+METAL[:2]+'_072219','grids/sigma_'+METAL[:2]+'_072219',setAsDefaults=True)
    for bv in [0.65]:#[1.5,1.75,1.9]:
        for li in np.linspace(-3.8,-5,5):#[1,2,3]:#[3.1,1]:
            my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=False,shadeScatter=False,errorbars=True,title='B-V= %s' % bv, bv_m=bv_m,upper_lim=upper_lim,metal_val=li)
            #baffles.baffles_age(bv,li=li,showPlots=True,pdfPage=pp)
            p = baf.get_posterior(bv,li,pdfPage=pp,showPlot=True,logPlot=False,upperLim = False,mamajekAge=True)
    pp.close()
    printName(name)
    

    """
    # Update default grids
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,rhk_fits,'grids/median_rhk_030719','grids/sigma_rhk_030719',True)
    
    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(bv_m,fits,upper_lim,'grids/median_li_071819','grids/sigma_li_071819',setAsDefaults=True)
    """


if  __name__ == "__main__":
    main()
