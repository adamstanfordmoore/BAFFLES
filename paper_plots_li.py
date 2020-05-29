"""
Adam Stanford-Moore
5/22/20
Code for plotting the lithium plots in the final BAFFLES paper
"""
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import fitting as my_fits
import probability as prob
import baffles
import li_constants as const
import plotting as my_plot
import readData
import fitting as my_fits
import utils
import bisect
from os.path import join
import os
if not os.path.exists('plots'):
    os.mkdir('plots')

METAL = "lithium"
bv_m, upper_lim, fits = readData.read_lithium()#fromFile=False,saveToFile=False)
bv_ca, ca_fits = readData.read_calcium()#fromFile=False,saveToFile=False)
print("Number of Li Stars= ",[len(x[0]) for x in bv_m])

def printName(n):
    print(" \n")

  
def main():
    #metal_vs_bv()
    #metal_vs_age()
    #metal_vs_age_subplots()
    #bldb()
    fit_hist()
    combined_validation_subplots()
    #moving_group()
    #notable_stars()
    #plot_fits_subplots()
    
    #---------Not in Paper-----------
    #combined_validation() 
    #get_CI_hyades_no_ul()
    #plot_fits()

    return


def moving_group():
    name = join('plots','moving_group_age.pdf')
    pp = PdfPages(name)
    baf = baffles.age_estimator('lithium')
    abdor = readData.abdor()
    fit = my_fits.poly_fit(abdor[0],abdor[1],2,scatter=True)
    plt.title("AB Dor")
    #my_plot.plot_fits([abdor],[fit],METAL,pdfPage=pp,showPlots=False)

    p = baf.posterior_product(abdor[0],abdor[1],None,abdor[2],pdfPage=pp,showPlot=False,
                            showStars=True,title='AB Dor',givenAge=149,givenErr=[-19,51])

    tuc = readData.tuchor()
    fit = my_fits.poly_fit(tuc[0],tuc[1],2,scatter=True)
    plt.title("Tuc/Hor")
    #my_plot.plot_fits([tuc],[fit],METAL,pdfPage=pp,showPlots=False)
    p2 = baf.posterior_product(tuc[0],tuc[1],None,tuc[2],pdfPage=pp,showPlot=False,
                                showStars=True,title='Tuc/Hor',givenAge=45,givenErr=4)

    print("We derive ages for AB Dor: $%d_{%d}^{+%d}$ Myr, and Tuc/Hor: $%d_{%d}^{+%d}$ Myr" \
           % (p.stats[2],p.stats[1]-p.stats[2],p.stats[3]-p.stats[2],
           p2.stats[2],p2.stats[1]-p2.stats[2],p2.stats[3]-p2.stats[2]))

    pp.close()
    printName(name)



def plot_fits(cluster=None):
    #Plot Fits
    name = join('plots', METAL + '_fits.pdf')
    pp=PdfPages(name)
    my_plot.plot_fits(bv_m,fits,METAL,pdfPage=pp,showPlots=False,upper_lim=upper_lim,
                    shadeScatter=True,specific_clusters=cluster)
    printName(name)
    pp.close()

def plot_fits_subplots():

    name = join('plots', METAL + '_fits_subplots.pdf')
    pp=PdfPages(name)
    fig,ax = plt.subplots(5,2,sharex=True,sharey=True,figsize=(12,12))

    plt.subplots_adjust(wspace=0, hspace=0)
    for i in range(len(bv_m)):
        color = const.COLORS[i]
        marker = const.MARKERS[i]
        
        my_plot.plot_fits_subplot(ax[int(i/2),i%2],i,bv_m,fits,METAL,upper_lim = upper_lim)
    
    # Set common labels
    fig.text(0.5, 0.04, r'$(B-V)_o$',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center')
    fig.text(0.06, 0.5, r'Log(Li EW/m$\AA$)',size=my_plot.AXIS_LABEL_SIZE, ha='center', 
            va='center', rotation='vertical')

    pp.savefig()
    plt.close()    
    printName(name)
    pp.close()



def metal_vs_bv():
    #Metal vs B-V
    name = join('plots',METAL + '_metal_vs_bv.pdf')
    pp=PdfPages(name)
    ###clusters = [0,4,6,7]
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=False,upper_lim=upper_lim,primordial_li=True,
                        specific_clusters = None,legend=False)
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=False,upper_lim=upper_lim,primordial_li=True,
                        fits_only=True,shadeScatter=False)
    printName(name)
    pp.close()

def metal_vs_age():
    #Metal vs Age
    name = join('plots', METAL + '_metal_v_age.pdf')
    pp=PdfPages(name)
    bvRange = const.BV_S 
    for bv in bvRange:
        print("bv: ",bv)
        my_plot.metal_vs_age(fits,METAL,bv,pp,showPlots=False,shadeScatter=False,\
                errorbars=False,title='B-V = %.3f' % bv, bv_m=bv_m,upper_lim=upper_lim,\
                logAge=True,plotStars=True)
    printName(name)
    pp.close()

def metal_vs_age_subplots():
    #Metal vs Age
    name = join('plots', METAL + '_metal_v_age_subplots.pdf')
    pp=PdfPages(name)
    fig,ax = plt.subplots(4,2,sharex=True,sharey=True,figsize=(10,12))
    #plt.tight_layout()
    bvRange = [.448,.621,.768,.965,1.137,1.310,1.506,1.703]
    for i,bv in enumerate(bvRange):
        print("bv: ",bv)
        my_plot.metal_vs_age_subplot(ax[int(i/2),i%2],fits,METAL,bv,pp,showPlots=False,shadeScatter=False,\
                errorbars=False,title='B-V = %.3f' % bv, bv_m=bv_m,upper_lim=upper_lim,\
                logAge=True,plotStars=True,legend=False)
    # Set common labels
    fig.text(0.5, 0.04, 'Age (Myr)',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center')
    fig.text(0.06, 0.5, r'Log(Li EW/m$\AA$)',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center', rotation='vertical')

    #plt.xlabel('Age (Myr)',size=my_plot.AXIS_LABEL_SIZE)
    #plt.ylabel(r'Log(Li EW/m$\AA$)',size=my_plot.AXIS_LABEL_SIZE)
    pp.savefig()
    plt.close()
    printName(name)
    pp.close()

def scatter_vs_bv():
    #Scatter vs B-V
    pp=PdfPages(join('plots', METAL + '_scatter_vs_bv_with_offset.pdf'))
    my_plot.scatter_vs_bv(fits,METAL,pp,showPlots=True)
    names.append(join('plots', METAL + '_scatter_vs_bv_with_offset.pdf'))
    pp.close()

def scatter_vs_age():
    #Scatter vs Age
    name = join('plots',METAL + '_scatter_vs_age.pdf')
    pp=PdfPages(name)
    my_plot.scatter_vs_age(fits,METAL,.65,pp,showPlots=True,bv_m=bv_m,upper_lim=upper_lim)
    pp.close()
    printName(name)

def fit_hist():
    #Fit histogram
    name = join('plots',METAL + '_hist_fit.pdf')
    pp=PdfPages(name)
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=False,upper_limits=upper_lim)
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=False,upper_limits=upper_lim,plot_cdf=True)
    printName(name)
    pp.close()
   
def bldb():
    #BLDB fit
    my_fits.bldb_fit(fits,plot=True)
    printName('bldb_vs_age.pdf')


def combined_validation():
    #Omitting each cluster
    name = join('plots',METAL + '_combined_validation.pdf')
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

def combined_validation_subplots():
    const = utils.init_constants(METAL)
    #Omitting each cluster
    name = join('plots',METAL + '_combined_validation_subplots.pdf')
    pp=PdfPages(name)
    baf_default = baffles.age_estimator(METAL)

    fig,ax = plt.subplots(3,2,figsize=(14,15))
    cmap = plt.cm.get_cmap('RdYlBu_r')
    norm = mpl.colors.Normalize(vmin=const.BV_RANGE[0], vmax=const.BV_RANGE[1])
    sc = plt.scatter([],[],c=[],norm=norm,cmap=cmap)
    #cbar = fig.colorbar(sc)

    fig.tight_layout(pad=.4,w_pad=1, h_pad=2)
    fig.subplots_adjust(left=0.06)
    fig.subplots_adjust(bottom=0.06)
    fig.subplots_adjust(top=.95)
    fig.subplots_adjust(right=0.9)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    fig.colorbar(sc, cax=cbar_ax)

    for index,i in enumerate([0,3,4,5,6,8]):    #[1,4,6,8]):
        print(const.CLUSTER_NAMES[i])
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p_val = baf.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=upper_lim[i])

        pl = ax[int(index/2),index%2]
        pl.plot(const.AGE,p_val.array,linewidth=2,linestyle='--',label='Posterior with removal')
        pl.set_title(const.CLUSTER_NAMES[i],size=my_plot.TITLE_SIZE)

        bv_arr = bv_m[i][0]
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=upper_lim[i],showStars=True)
        #print(p.stars_posteriors)
        age,y = const.AGE,p.array
        givenAge=const.CLUSTER_AGES[i]

        for star,post in enumerate(p.stars_posteriors):
            color = cmap(norm(bv_arr[star]))
            prob.scale_to_height(post,np.max(y))
            pl.plot(const.AGE,post,alpha = 1,linewidth=1,color=color,zorder=0)

        if (givenAge):
            print('Isochronal age exists within %f %% CI' % prob.get_percentile(age,y,givenAge))
            pl.axvline(x=givenAge,color='r',label='Isochronal age: %d Myr' % givenAge)
            #if (givenErr):
            #    if (type(givenErr) == float or type(givenErr) == int):
            #        givenErr = [-1*givenErr,givenErr]
            #    pl.axvspan(givenAge+givenErr[0], givenAge+givenErr[1], alpha=0.2, color='r',zorder=0)
        #if (mamajekAge):
        #    plt.axvline(x=mamajekAge,color='C2',label='MH08 age: %d' % mamajekAge)


        pl.plot(age,y,color = 'C0',linewidth=2)
        stat = p.stats
        age2 = np.linspace(stat[0],stat[-1],500)
        interp = my_fits.piecewise(age,y)
        y2 = interp(age2)
        pl.vlines(x=stat[2],ymin= 0,ymax= interp(stat[2]), \
                label='BAFFLES median age: %.3g Myr' % stat[2] ,color = 'orange')
        pl.fill_between(age2,y2, where= (age2 >= stat[1]) & (age2 <= stat[3]),color='.3', \
                label='68%% CI: %.2g - %.2g' % (stat[1],stat[-2]))
        pl.fill_between(age2,y2, where= (age2 >= stat[0]) & (age2 <= stat[-1]),color='.6',\
                alpha=0.5, label='95%% CI: %.2g - %.2g' % (stat[0],stat[-1]))
        pl.set_ylim([0,np.max(y)*1.5])
        r = my_plot.getAgeRange(p.stats,p.stars_posteriors,givenAge)
        pl.set_xlim(r)
        pl.legend()
        pl.minorticks_on()
        pl.tick_params(axis='both',which='both',right=True,top=True)

    # Set common labels
    fig.text(0.5, 0.02, 'Age (Myr)',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center')
    fig.text(0.01, 0.5, 'Probability Density (Myr^-1)',size=my_plot.AXIS_LABEL_SIZE, ha='center', 
            va='center', rotation='vertical')
    fig.text(0.99, 0.5, 'B-V',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center', rotation='vertical')
    pp.savefig()
    plt.close()
    printName(name)
    pp.close()

def get_CI():
    const = utils.init_constants(METAL)
    baf_default = baffles.age_estimator(METAL)

    for index,i in enumerate([1,2,7,9]):  
        print("\n", const.CLUSTER_NAMES[i])

        bv_arr = bv_m[i][0]
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=upper_lim[i],showStars=True)
        age,y = const.AGE,p.array
        givenAge=const.CLUSTER_AGES[i]

        print('Isochronal age exists within %f %% CI' % prob.get_percentile(age,y,givenAge))

def get_CI_hyades_no_ul():
    const = utils.init_constants(METAL)
    baf_default = baffles.age_estimator(METAL)
 
    mask =  (~np.array(upper_lim[8])) | (bv_m[8][0] > 0.55) 
    print("Total stars:",len(bv_m[8][0]))
    print("Stars selsected=:",np.sum(mask))
    b,l,ul = bv_m[8][0][mask],bv_m[8][1][mask],np.array(upper_lim[8])[mask]
    baf_default.posterior_product(b,l,upperLim_arr=ul,showPlot=True,showStars=True,title='Hyades B-V > .55',givenAge=700)
    
    #bv_arr = bv_m[i][0][~np.array(upper_lim[i])]
    #metal_array = bv_m[i][1][~np.array(upper_lim[i])]
    #ul_arr = upper_lim[i][~np.array(upper_lim[i])] #np.array(upper_lim[i])[np.array(bv_m[i][0]) > 0.55]
    #p = baf_default.posterior_product(bv_arr,metal_array,upperLim_arr=ul_arr,showPlot=True,title='Hyades B-V > .55',givenAge=700)



def omitting(validation=False):
    #Omitting each cluster
    #name = join('plots',METAL + '_self_validation.pdf')
    name = join('plots',METAL + '_omit_clusters.pdf') if not validation else \
           join('plots',METAL + '_self_validation.pdf')
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
    name = join('plots','notable_stars.pdf')
    names = ["HR 2562","HD 206893","TW PsA"]
    bv = [.45,.44,1.1]
    bv_err = [np.sqrt(.014**2+.01**2),np.sqrt(.02**2+.01**2),.03]   # .03 b-v error for TW PsA?
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

    my_plot.metal_vs_bv(bv_ca,ca_fits,'calcium',None,False,specific_clusters=[0,5,7,8],
                        legend=False,textlabels=True)
    plt.plot(bv,rhk,marker='s',markersize=markerSize,color='w',linestyle='None',zorder=9)
    for i in [0,1]:
        plt.plot(bv[i],rhk[i],marker=markers[i],markersize=markerSize,color=colors[i],
                linestyle='None',zorder=10,label=names[i])
    plt.legend()
    plt.xlim([.42,.9])
    pp.savefig()
    #plt.show()
    plt.close()

    my_plot.metal_vs_bv(bv_m,fits,'lithium',None,False,upper_lim=upper_lim,specific_clusters=[0,1,4,6,8,9])
    plt.plot(bv,li,marker='s',markersize=markerSize,color='w',linestyle='None',zorder=9)
    for i in [0,1,2]:
        plt.plot(bv[i],np.log10(li[i]),marker=markers[i],markersize=markerSize,color=colors[i],
                linestyle='None',zorder=10,label=names[i])
    plt.legend()
    pp.savefig()
    #plt.show()
    plt.close()

    baf_ca = baffles.age_estimator('calcium')
    baf_li = baffles.age_estimator('lithium')
    for i in [0,1]:
        print(names[i])
        plt.plot([],[],'C0',linewidth=2,label='Final Age Posterior')
    
        p_li = baf_li.get_posterior(bv[i],li[i],bv_uncertainty=bv_err[i],measure_err=li_err[i],upperLim=False)
        p_ca = baf_ca.get_posterior(None,rhk[i])
        product = prob.normalize(const.AGE,p_ca.array*p_li.array)
        prod_stats=prob.stats(const.AGE,product)
        my_plot.posterior(const.AGE, product, prod_stats,names[i],logPlot=False)
        plt.plot(const.AGE, p_ca.array,color='C3',label="Calcium Posterior")
        plt.plot(const.AGE, p_li.array,color='C2',label="Lithium Posterior")
        plt.axvspan(age_range[i][0],age_range[i][1], alpha=0.2, color='r',
                    label=r'Literature age: %d - %d Myr' % tuple(age_range[i]),zorder=0)
        plt.axvline(x=mamaAge[i],color='C5',linestyle='--',label='MH08 age: %d' % mamaAge[i])
        #plt.xlim([0,490])
        plt.legend()
        pp.savefig()
        #plt.show()
        plt.close()
        print("%d Myr (68\\%%CI: %d - %d Myr)" % (utils.round_sigs(p_ca.stats[2],2),
              utils.round_sigs(p_ca.stats[1],2),utils.round_sigs(p_ca.stats[3],2)))
        print("%.1f Gyr (68\\%%CI: %.1f - %.1f Gyr)" % (p_li.stats[2]/1000,p_li.stats[1]/1000,p_li.stats[3]/1000))
        print("%d Myr, with a 68\\%% confidence interval between %d Myr - %d Myr" % (utils.round_sigs(prod_stats[2],2),
               utils.round_sigs(prod_stats[1],2),utils.round_sigs(prod_stats[3],2)))



    print("TW PsA")
    plt.axvspan(age_range[-1][0],age_range[-1][1], alpha=0.2, color='r',zorder = 0)
    plt.axvspan(360-140,360+140, alpha=0.2, color='b',zorder=0)
    p_li = baf_li.get_posterior(bv[2],li[2],bv_uncertainty=bv_err[2],measure_err=li_err[2],upperLim=False)
    print("we report an age of %d Myr with a 68\\%% confidence interval between %d Myr - %d Myr\
          (third panel of Fig. \\ref{fig:notable_stars}), consistent with Mamajek's lithium age,\
           but a factor of $\sim$%.1f too young for his final adopted age." % (p_li.stats[2],p_li.stats[1],p_li.stats[3],
          440/p_li.stats[2]))
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
    print("to get a final age for the system, $%d^{+%d}_{%d}$ Myr." % (stat[2],stat[3]-stat[2],stat[1]-stat[2]))
    plt.vlines(x=stat[2],ymin= 0,ymax= y[bisect.bisect_left(const.AGE,stat[2])], \
                    label='Final median age: %.3g Myr' % stat[2] ,color = 'orange')
    plt.fill_between(const.AGE,y, where= (const.AGE >= stat[1]) & (const.AGE <= stat[3]),color='.4', \
                        label='68%% CI: %.2g - %.2g' % (stat[1],stat[-2]))
    plt.title(names[2],size=my_plot.TITLE_SIZE)
    plt.xlim([0,1200])
    plt.legend()
    plt.ylabel('Probability Density (Myr^-1)',size=my_plot.AXIS_LABEL_SIZE)
    plt.xlabel("Age (Myr)",size=my_plot.AXIS_LABEL_SIZE)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    pp.savefig()
    #plt.show()
    plt.close()

    printName(name)
    pp.close()

def robs_fomalhaut():
    import astropy
    from scipy.signal import savgol_filter
    t = astropy.io.fits.open(join('data','Fom-age-pdf.fits'))
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


if  __name__ == "__main__":
    main()
