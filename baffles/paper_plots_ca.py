"""
Adam Stanford-Moore
5/22/20
Code for plotting the calcium plots in the final BAFFLES paper
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
import ca_constants as const
import plotting as my_plot
import readData
import utils
from os.path import join
import os
if not os.path.exists('plots'):
    os.mkdir('plots')

METAL = "calcium"
upper_lim = None
bv_m,fits = readData.read_calcium()#fromFile=False,saveToFile=False,fit_degree=0)
print("Num Calcium Stars= ", [len(x[0]) for x in bv_m])

def main(): 
    #metal_vs_bv()
    #metal_vs_age()
    #scatter_vs_age()
    #fit_hist()
    #baffles_vs_mamajek()
    #combined_validation_subplots()
    

    #---------Not in Paper-----------
    #plot_fits()
    #combined_validation()
    #posteriors()
    #nearest_stars_hist()
    #print_BIC()
    return


def print_BIC():
    bv_m,fits = readData.read_calcium(fromFile=False,saveToFile=False,fit_degree=0)
    dof = len(bv_m) #1 parameter for each cluster
    BIC1 = my_fits.get_fit_BIC(bv_m,fits,dof)
    print("BIC constant fits = ", BIC1)
    

    dof = 2*len(bv_m) # 2 for each cluster
    bv_m,fits = readData.read_calcium(fromFile=False,saveToFile=False,fit_degree=1)
    BIC2 = my_fits.get_fit_BIC(bv_m,fits,dof)
    print("BIC linear fits = ", BIC2)
    #for i in range(len(bv_m)):
    #    plt.scatter(bv_m[i][0],bv_m[i][1])
    #    plt.plot(bv_m[i][0],fits[i][0](bv_m[i][0]))
    #    plt.show()
    print("Delta BIC = ", BIC2 - BIC1)

def printName(n):
    print("-----\n")


def plot_fits(cluster_indices=None):
    #Plot Fits
    name = join('plots', METAL + '_fits.pdf')
    pp=PdfPages(name)
    my_plot.plot_fits(bv_m,fits,METAL,pdfPage=pp,showPlots=False,upper_lim=upper_lim,
                        specific_clusters=cluster_indices)
    pp.close()
    printName(name)

def metal_vs_bv():
    #Metal vs B-V
    name = join('plots',METAL + '_metal_vs_bv.pdf')
    pp=PdfPages(name)
    ##clusters = [0,4,6,7]
    my_plot.metal_vs_bv(bv_m,fits,METAL,pp,showPlots=False,legend=False)#,specific_clusters = clusters)
    
    bv_m_linear,fits_linear = readData.read_calcium(fromFile=False,saveToFile=False,fit_degree=1)
    my_plot.metal_vs_bv(bv_m_linear,fits_linear,METAL,pp,showPlots=False)#,specific_clusters = clusters)
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
                errorbars=True,title=' ', bv_m=bv_m,upper_lim=upper_lim,\
                logAge=True,plotStars=False,mamajek_poly=True)
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
    my_plot.scatter_vs_age(fits,METAL,.65,pp,showPlots=False,bv_m=bv_m,upper_lim=upper_lim)
    pp.close()
    printName(name)

def fit_hist():
    #Fit histogram
    name = join('plots',METAL + '_hist_fit.pdf')
    pp=PdfPages(name)
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=False)
    my_plot.fit_histogram(bv_m,fits,METAL,pp,showPlots=False,plot_cdf=True)
    printName(name)
    pp.close()

def baffles_vs_mamajek():
    # BAFFLES vs Mamajek
    name = join('plots','baffles_vs_mamajek.pdf')
    pp = PdfPages(name)
    for i in range(len(bv_m)):
        my_plot.baffles_vs_mamajek(bv_m,fits,i,pdfPage=pp,showPlots=False,title=None,mamaProduct=False)
    pp.close()
    printName(name)

def omitting(validation=False):
    #Omitting each cluster
    #name = join('plots',METAL + '_self_validation.pdf')
    name = join('plots',METAL + '_omit_clusters.pdf') if not validation else \
            join('plots',METAL + '_self_validation.pdf')
    pp = PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=validation)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        #mamaProductAge = utils.getMamaProductAge(bv_m[i][1])
        #print(mamaProductAge
        if not validation: baf.make_grids(bv_m,fits,omit_cluster=i)
        p = baf.posterior_product(bv_m[i][0],bv_m[i][1],pdfPage=pp,showPlot=False,\
                showStars=True,givenAge=const.CLUSTER_AGES[i],\
                title= const.CLUSTER_NAMES[i])
        #print("True: %.3g BAFFLES: %.3g Mama: %.3g" % (const.CLUSTER_AGES[i],p.stats[2],mamaProductAge)
    printName(name)
    pp.close()


def combined_validation():
    #Omitting each cluster
    name = join('plots',METAL + '_combined_validation.pdf')
    pp=PdfPages(name)
    baf_default = baffles.age_estimator(METAL)
    for i in range(len(bv_m)):
        print(const.CLUSTER_NAMES[i])
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p_val = baf.posterior_product(bv_m[i][0],bv_m[i][1])

        plt.plot(const.AGE,p_val.array,linewidth=2,linestyle='--',label='Posterior with removal')
        #Since given age provided, prints where isochronal age lies without removal
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],\
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

    fig.tight_layout(pad=.4,w_pad=1, h_pad=2)
    fig.subplots_adjust(left=0.06)
    fig.subplots_adjust(bottom=0.06)
    fig.subplots_adjust(top=.95)
    fig.subplots_adjust(right=0.9)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    fig.colorbar(sc, cax=cbar_ax)

    for index,i in enumerate([1,3,5,6,7,8]):        #[0,4,5,6,7,8] for submission 2!
        print(const.CLUSTER_NAMES[i])
        baf = baffles.age_estimator(METAL,default_grids=False)
        baf.make_grids(bv_m,fits,omit_cluster=i)
        p_val = baf.posterior_product(bv_m[i][0],bv_m[i][1])

        pl = ax[int(index/2),index%2]
        pl.plot(const.AGE,p_val.array,linewidth=2,linestyle='--',label='Posterior with removal')
        pl.set_title(const.CLUSTER_NAMES[i],size=my_plot.TITLE_SIZE)

        bv_arr = bv_m[i][0]
        p = baf_default.posterior_product(bv_m[i][0],bv_m[i][1],upperLim_arr=None,showStars=True)
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
    fig.text(0.01, 0.5, 'Probability Density (Myr^-1)',size=my_plot.AXIS_LABEL_SIZE, 
            ha='center', va='center', rotation='vertical')
    fig.text(0.99, 0.5, 'B-V',size=my_plot.AXIS_LABEL_SIZE, ha='center', va='center', 
            rotation='vertical')
    pp.savefig()
    plt.close()
    printName(name)
    pp.close()


def posteriors():
    #Making posteriors
    name = join('plots', METAL + '_posteriors.pdf')
    pp = PdfPages(name)
    baf = baffles.age_estimator(METAL,default_grids=False)
    baf.make_grids(bv_m,fits,upper_lim)#,omit_cluster=0)
    for bv in [0.65]:
        for li in np.linspace(-3.8,-5,5):
            p = baf.get_posterior(bv,li,pdfPage=pp,showPlot=True,logPlot=False,upperLim = False,
                                  mamajekAge=True)
    printName(name)
    pp.close()


def nearest_stars_hist():
    t = np.genfromtxt(join("data","MH08_table13.txt"),delimiter=';',dtype=str,skip_header=75)
    bv,rhk = t[:,6].astype(np.float),t[:,8].astype(np.float)
    #mask = (const.BV_RANGE[0] <= bv) & (bv <= const.BV_RANGE[1]) & \
    mask = (.8 <= bv) & (bv <= .9) & \
           (const.METAL_RANGE[0] <= rhk) & (rhk <= const.METAL_RANGE[1])
    bv,rhk = bv[mask],rhk[mask]
    baf = baffles.age_estimator('calcium')
    for i in range(1):
        ages = []#,array = [],np.zeros(1000)
        for b,r in zip(bv,rhk):
            p = baf.get_posterior(b,r)
            #chance = p.array/np.sum(p.array)
            a = p.stats[2]#np.random.choice(const.AGE,p=chance)
            ages.append(a)
            #array += p.array
            #plt.hist(ages)
            #plt.xlabel("Age (Myr)")
            #plt.ylabel("Frequency")
            #plt.show()
        plt.hist(ages)
        plt.xlabel("Age (Myr)")
        plt.ylabel("Frequency")
        plt.title(".8 <= B-V <= .9")
        plt.show()

        #x,cdf = prob.hist_cdf(ages)
        #plt.plot(x,cdf,color='gray')

    #plt.show()

        #pp.close()
        #prob.normalize(const.AGE,array)
        #plt.plot(const.AGE,array)
        #plt.show()


if  __name__ == "__main__":
    main()
