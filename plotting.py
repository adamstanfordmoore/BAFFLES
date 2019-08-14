import matplotlib.pyplot as plt
import bisect
import ca_constants as const
import numpy as np
import probability as prob
import fitting as my_fits
from scipy import interpolate
import baffles
import utils
from scipy.stats import lognorm,norm
#COLORS = ['C0','C4','C2','C1','C5','C6','C7','C8']

import matplotlib as mpl
import matplotlib.cm as cm

#givenErr = [-5,5] or 5
def posterior(age,y,stat,title=' ',pp=None,showPlot=False,starArray = [],\
        givenAge=None,givenErr = None,mamajekAge=None,logPlot=False,bv_arr = None,\
        metal='calcium'):
    const = init_constants(metal)
    isUpperLim = len(stat) == 4 #arbitrarily [.002,.05,.34,1]
    plt.title(title)
    
    cmap,norm,sc = None,None,None
    if bv_arr is not None:
        cmap = plt.cm.get_cmap('RdYlBu_r')
        #cmap = plt.cm.get_cmap('bwr')
        norm = mpl.colors.Normalize(vmin=const.BV_RANGE[0], vmax=const.BV_RANGE[1])#np.max(bv_arr))
    for i,post in enumerate(starArray):
        color = '.5'
        if bv_arr is not None: 
            color = cmap(norm(bv_arr[i]))
        
        prob.scale_to_height(post,np.max(y))
        plt.plot(const.AGE,post,alpha = 1,linewidth=1,color=color,zorder=0)
    
    if bv_arr is not None:
        sc = plt.scatter([],[],c=[],norm=norm,cmap=cmap)
        cbar = plt.colorbar(sc)
        cbar.ax.set_ylabel('B-V', size = 18,rotation=270,labelpad=17)

    if (givenAge):
        plt.axvline(x=givenAge,color='r',label='Isochronal age: %d Myr' % givenAge)
        if (givenErr):
            if (type(givenErr) == float or type(givenErr) == int):
                givenErr = [-1*givenErr,givenErr]
            plt.axvspan(givenAge+givenErr[0], givenAge+givenErr[1], alpha=0.2, color='r',zorder=0)
    if (mamajekAge):
        plt.axvline(x=mamajekAge,color='C2',label='MH08 age: %d' % mamajekAge)
    if (logPlot or isUpperLim):
        plt.semilogx(age,y,color = 'C0',linewidth=2)
    else:
        plt.plot(age,y,color = 'C0',linewidth=2)
    shadeStats(age,y,stat,isUpperLim)
    plt.xlabel(u'Age (Myr)',size=18)
    plt.ylabel('Probability',size=18)
    plt.legend() #loc="upper right")
    if (not logPlot and not isUpperLim):
        r = getAgeRange(stat,starArray)
        plt.xlim(r)
    if (len(starArray) > 1):
            plt.ylim([0,np.max(y)*1.5])
    plt.tight_layout()
    if (pp):
        pp.savefig()
    if (showPlot):
        plt.show()
    if pp is not None:
        plt.close()

def shadeStats(age,y,stat,upperLim):
    if upperLim:
        plt.fill_between(age,y, where= (age >= stat[2]) & (age <= stat[3]),\
                color='.3', label=r'1$\sigma$ lower lim: %.3g Myr' % stat[2])
        plt.fill_between(age,y, where= (age >= stat[1]) & (age <= stat[2]),\
                color='.6',alpha=1, label=r'2$\sigma$ lower lim: %.3g Myr' % stat[1])
        plt.fill_between(age,y, where= (age >= stat[0]) & (age <= stat[1]),\
                color='.9',alpha=1,label=r'3$\sigma$ lower lim: %.3g Myr' % stat[0])
        return
    
    age2 = np.linspace(stat[0],stat[-1],500)
    interp = my_fits.piecewise(age,y)
    y2 = interp(age2)

    plt.vlines(x=stat[2],ymin= 0,ymax= interp(stat[2]), \
            label='BAFFLES median age: %.3g Myr' % stat[2] ,color = 'orange')
    plt.fill_between(age2,y2, where= (age2 >= stat[1]) & (age2 <= stat[3]),color='.3', \
            label='68%% CI: %.2g - %.2g' % (stat[1],stat[-2]))
    plt.fill_between(age2,y2, where= (age2 >= stat[0]) & (age2 <= stat[-1]),color='.6',\
            alpha=0.5, label='95%% CI: %.2g - %.2g' % (stat[0],stat[-1]))
    #plt.vlines(x=stat[2],ymin= 0,ymax= y[bisect.bisect_left(const.AGE,stat[2])], \
    #        label='BAFFLES Age: %.3g Myr' % stat[2] ,color = 'orange')
    #plt.fill_between(age,y, where= (age >= stat[1]) & (age <= stat[3]),color='.3', \
    #        label='68%% CI: %.2g - %.2g' % (stat[1],stat[-2]))
    #plt.fill_between(age,y, where= (age >= stat[0]) & (age <= stat[-1]),color='.6',\
    #        alpha=0.5, label='95%% CI: %.2g - %.2g' % (stat[0],stat[-1]))

#Determines the domain for plotting posterior in linear space
def getAgeRange(stat,starArray):
        l = 0
        u = const.GALAXY_AGE
        median = stat[2]
        sigma = (stat[3] - stat[1])/2
        #sigma = median - stat[1]
        #sigma = stat[3] - median
        if (len(starArray) > 1):
            if (sigma < 10):
                sigma = 10
            #plt.xlim([r[0] - 2*sig,2*sig + r[1]])
            #plt.xlim([0,80])
        if (median - 4*sigma > l):
            l = median - 4*sigma
        if (median + 4*sigma < u):
            u = median + 4*sigma
        return [l,u]


def shade_scatter(fit,BV):
    y_top,y_bottom = fit[0](BV), fit[0](BV)
    y_top += fit[1](BV)
    y_bottom -= fit[1](BV)
    plt.fill_between(BV,y_bottom,y_top,alpha=.3)

def plot_fits(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,shadeScatter=True,specific_clusters=None):
    if (not pdfPage and not showPlots):
        return
    const = init_constants(metal)   
    for i in range(len(fits)):
        if specific_clusters is not None and i not in specific_clusters: continue
        ax = plt.gca()
        color = next(ax._get_lines.prop_cycler)['color']
        for j in range(len(bv_m[i][0])):
            if (upper_lim is not None and upper_lim[i][j]):
                plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=color,marker=const.DOWN_ARROW)
            else:
                plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=color)

        #l = 'piecewise: %d' % SEG if converged[i][s] else 'Not Converged: %d' % SEG
        #l = 'piecewise, bin size: %d' % SEG if converged[i][s] else 'Not Converged, bin size: %d' % SEG
        plt.plot(const.BV,fits[i][0](const.BV),dashes=[2,2], color=color)
        
        if shadeScatter:
            shade_scatter(fits[i],const.BV.tolist())
        plt.title(const.CLUSTER_NAMES[i])
        plt.axis(const.BV_RANGE + const.METAL_RANGE)
        plt.xlabel(r'$(B-V)_o$',size=18)
        set_ylabel(metal)
        #plt.legend()
        plt.tight_layout()
        if (pdfPage):
            pdfPage.savefig()
        if (showPlots):
            plt.show()
        plt.close()
        """
        resid = residuals(bv_m[i][0],bv_m[i][1],fits[i][s][0])
        plt.scatter(bv_m[i][0],resid,label = fit_names[s])
        plt.axhline(y=0,linestyle='--',color='k')
        plt.title(CLUSTER_NAMES[i] + ' Residuals for LogEW vs. ' + r'$(B-V)_o$')
        plt.xlim(BV_RANGE)
        plt.xlabel(r'$(B-V)_o$',size=18)
        plt.ylabel(u'logEW',size=18)
        plt.legend()
        pp.savefig()
        #plt.show()
        plt.close()
        """

def metal_vs_bv(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,shadeScatter=False,title=None,primordial_li = False,fits_only=False,specific_clusters = None,legend=True):
    const = init_constants(metal)
    plt.xlabel(r'$(B-V)_0$',size=18)
    #plt.ylabel(r'EW Li (m$\AA$)',size=18)
    set_ylabel(metal)
    
    plt.axis([const.BV_RANGE[0]-.01,const.BV_RANGE[1]] + const.METAL_RANGE)
    #plt.axis(const.BV_RANGE + const.METAL_RANGE)
    #plt.axis(const.BV_RANGE + np.power(10,const.METAL_RANGE).tolist())
    #plt.axis([const.BV_RANGE[0],1, 0,400])
    
    #plot primordial lithium
    if (primordial_li and metal.lower()[0]=='l'):
        #ngc2264_fit = fits[const.CLUSTER_NAMES.index("NGC2264")][0]
        primordial_li_fit = my_fits.MIST_primordial_li()#ngc2264_fit,fromFile=False,saveToFile=True) 
        plt.plot(const.BV,primordial_li_fit(const.BV),color=const.PRIM_COLOR,label="Primordial EW Li")
        #plt.plot(const.BV,primoIrdial_li_fit(const.BV),color=const.PRIM_COLOR,label="Primordial EW Li")
    
    clusters_to_plot = specific_clusters if specific_clusters != None else range(len(bv_m))
    for c in clusters_to_plot:
        ax = plt.gca()
        #color = const.COLORS[c] #next(ax._get_lines.prop_cycler)['color']
        if (not fits_only):
            for i in range(len(bv_m[c][0])):
                if (upper_lim and upper_lim[c][i]):
                    plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=const.COLORS[c],marker=const.DOWN_ARROW)
                else:
                    plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=const.COLORS[c],marker=const.MARKERS[c])
        if legend:
            plt.scatter([],[],color=const.COLORS[c],marker=const.MARKERS[c],label=const.CLUSTER_NAMES[c])
        
        if fits_only:
            plt.plot(const.BV,fits[c][0](const.BV), color=const.COLORS[c])
        else:
            plt.plot(const.BV,fits[c][0](const.BV),dashes=[2,2], color=const.COLORS[c])
        
        if (shadeScatter):
            shade_scatter(fits[c],const.BV.tolist())

    if not legend and metal == 'calcium':
        plt.text(.7,-4.005,"Sco-Cen",size=13,color=const.COLORS[0])
        plt.text(.725,-4.27,"Pleiades",size=13,color=const.COLORS[4])
        plt.text(.75,-4.55,"Hyades",size=13,color=const.COLORS[6])
        plt.text(.61,-4.77,"M67",size=13,color=const.COLORS[7])

    if (title):
        plt.title(title,size=18)
        #else:
        #   plt.title(metal[0].upper() + metal[1:] + r' Fits')
    
    if legend:
        leg = plt.legend(loc='lower right')
        for i in range(len(leg.get_texts())):
            if (primordial_li and metal == 'lithium'):
                if (i==0):
                    plt.setp(leg.get_texts()[0],color=const.PRIM_COLOR)
                else:
                    plt.setp(leg.get_texts()[i],color=const.COLORS[clusters_to_plot[i-1]])
            else:
                plt.setp(leg.get_texts()[i],color=const.COLORS[clusters_to_plot[i]])
        
    plt.tight_layout()
        
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()

#RHK vs age polynomial
def metal_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,\
        shadeScatter=False,errorbars=True,bv_m=None,upper_lim=None, \
        mamajek_poly=False,metal_val=None,logAge = True,plotStars=False,omit_cluster=None):
    const = init_constants(metal)
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = my_fits.get_valid_metal(bv,fits,const,None,None,omit_cluster)

    

    if (metal[0].lower()=='l'):
        plt.figure(figsize=(7,5))
    if metal_val:
        plt.hlines(metal_val,0,const.GALAXY_AGE,linestyle='dashed',color='orange',\
                label=metal[:2].title() + ' = %.2f' % metal_val)
    
    
    MARKERS,COLORS = const.MARKERS,const.COLORS
    if metal=='lithium':
        if bv >= const.BLDB_LOWER_LIM: 
            colors = [const.COLORS[const.CLUSTER_NAMES.index(name)] for name in CLUSTER_NAMES[1:-1]]
            markers = [const.MARKERS[const.CLUSTER_NAMES.index(name)] for name in CLUSTER_NAMES[1:-1]]    
            MARKERS = [const.PRIM_MARKER] + markers + [const.BLDB_MARKER]
            COLORS = [const.PRIM_COLOR] + colors + [const.BLDB_COLOR]
        else:
            colors = [const.COLORS[const.CLUSTER_NAMES.index(name)] for name in CLUSTER_NAMES[1:]]
            markers = [const.MARKERS[const.CLUSTER_NAMES.index(name)] for name in CLUSTER_NAMES[1:]]    
            MARKERS = [const.PRIM_MARKER] + markers
            COLORS = [const.PRIM_COLOR] + colors

    for i in range(len(rhk)):
        if (errorbars):
            plt.errorbar(CLUSTER_AGES[i], rhk[i], yerr=scatter[i],color=COLORS[i],\
                    marker=MARKERS[i],markersize=8,capsize=4,zorder=10)
            plt.scatter([], [], color=COLORS[i],marker=MARKERS[i],label=CLUSTER_NAMES[i])
        else:
            plt.scatter(CLUSTER_AGES[i], rhk[i],color=COLORS[i],marker=MARKERS[i],\
                    label=CLUSTER_NAMES[i],s=80,zorder=10)
    
    li_scatter_fit = my_fits.fit_two_scatters(bv_m,fits,upper_lim) if metal == 'lithium' else None
    mu,sig,lbl = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,metal,li_scatter_fit,omit_cluster)

    plt.plot(const.AGE,mu(const.AGE),label=lbl)
    
    if (plotStars and bv_m is not None):
        bv_threshold = 0.05 if metal=='lithium' else 1
        for c in range(len(bv_m)):
            bv_arr = np.array(bv_m[c][0])
            mask = (bv_arr >= bv - bv_threshold) & (bv_arr <= bv + bv_threshold)
            li = np.array(bv_m[c][1])[mask]
            star_ages = [const.CLUSTER_AGES[c]]*len(li)
            ul = np.array(upper_lim[c])[mask] if upper_lim is not None else None
            for j in range(len(li)):
                m = const.DOWN_ARROW if (upper_lim is not None and ul[j]) else '+'
                color = const.COLORS[c]# if metal == 'lithium' else const.COLORS[c]
                plt.scatter(const.CLUSTER_AGES[c],li[j],marker=m,s=30,color=color) #for primordial li

    if (shadeScatter):
        plt.fill_between(const.AGE,mu(const.AGE) - sig(const.AGE),\
                mu(const.AGE) + sig(const.AGE),alpha=.2,color='C0')
        #plt.fill_between(const.AGE,fit(np.log10(const.AGE)) - sig,fit(np.log10(const.AGE)) + sig,alpha=.2,color='C0')
        #plt.scatter(CLUSTER_AGES[i], rhk[i],color=const.COLORS[i],marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    
    if (mamajek_poly and metal[0].lower()=='c'):
        plt.plot(const.AGE,utils.getMamaRHK(const.AGE),linestyle='--',color='gray',\
                label='M & H 2008')
    
    ax = plt.gca()
    if logAge: ax.set_xscale('log')
    plt.axis([1,const.GALAXY_AGE]+const.METAL_RANGE)
    plt.legend()
    plt.xlabel('Age (Myr)',size=18)
    ylabel = set_ylabel(metal)
    if (title):
        plt.title(title,size=18)
    #else:
        #   plt.title('Mean ' + ylabel +  ' per Cluster at $(B-V)_o$ = %.2f' % bv,size=18)
    plt.tight_layout()
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    plt.close()        

def scatter_vs_bv(fits,metal,pdfPage=None,showPlots=False,title=None):
    const = init_constants(metal)
    plt.xlabel(r'$(B-V)_0$',size=18)
    #plt.ylabel(r'EW Li (m$\AA$)',size=18)
    m = 'Log(EWLi/m$\AA$)'
    plt.ylabel('Scatter in ' + m,size=16)
    #plt.axis(const.BV_RANGE + np.power(10,const.METAL_RANGE).tolist())
    #plt.axis([const.BV_RANGE[0],1, 0,400])
    for c in range(len(fits)):
        ax = plt.gca()
        color = const.COLORS[c] #next(ax._get_lines.prop_cycler)['color']
        plt.plot(const.BV,fits[c][1](const.BV),dashes=[2,2],label=const.CLUSTER_NAMES[c], color=color)
    plt.legend()
    if (title):
        plt.title(title,size=18)
        #else:
        #   plt.title(metal[0].upper() + metal[1:] + r' Fits')
    plt.tight_layout()
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    #plt.close()


def scatter_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,bv_m=None,upper_lim=None):
    const = init_constants(metal)
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = my_fits.get_valid_metal(bv,fits,const)
    
    for i in range(len(scatter)):
        plt.scatter(CLUSTER_AGES[i], scatter[i],color=const.COLORS[i],\
                marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    
    li_scatter_fit = my_fits.fit_two_scatters(bv_m,fits,upper_lim) if metal=='lithium' else None
    mu,sig,lbl = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,metal,li_scatter_fit)

    #plt.plot(const.AGE,mu(const.AGE),label=lbl)
   
    lbl = 'gaussian fit' if metal=='calcium' else 'meas + astro' #r'constant $\sigma$ = %.3f' % fit(.65) 
    plt.plot(const.AGE,sig(const.AGE),'--',label=lbl,color='orange')

    ax = plt.gca()
    ax.set_xscale('log')
    plt.legend()
    plt.xlabel('Age (Myr)',size=18)
    x_axis  = 'Log(EWLi/m$\AA$)'
    if (metal.lower()[0] == 'c'):
        x_axis = "Log(R'" + r'$_{HK})$'
    plt.ylabel('Scatter in ' + x_axis ,size=16)
    if (title):
        plt.title(title,size=18)
    #else:
        #   plt.title('Mean ' + ylabel +  ' per Cluster at $(B-V)_o$ = %.2f' % bv,size=18)
    plt.tight_layout()
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    plt.close()        


def set_ylabel(metal):
    if (metal[0].lower() == 'c'):
        #plt.ylabel(u'logR\'HK',size=18)
        plt.ylabel("Log(R'" + r'$_{HK})$',size=16)
        #plt.ylabel(u'Log(R\'HK) -- Ca Emission Strength',size=16)
        return 'Calcium Emission Strength'
    elif (metal[0].lower() == 'l'):
        #plt.ylabel('Lithium Abundance',size=16)
        #plt.ylabel(u'log(EW/mA)',size=18)
        plt.ylabel(r'Log(EWLi/m$\AA$)',size=18)
        return 'Lithium Abundance'

def init_constants(metal):
    if (metal[0].lower() == 'c'):
        import ca_constants as const
    elif (metal[0].lower() == 'l'):
        import li_constants as const
    else:
        raise RuntimeError("No metal specified. Please enter lithium or calcium")
    return const

"""
#plots gaussian cluster scatter at a given bv value
def plot_cluster_scatter(bv,rhk_fits,scatter_fits,RHK):
        plt.ylabel(r'Probability',size=18)
        plt.xlabel(u'logR\'HK',size=18)
        for i in range(len(CLUSTER_NAMES)):
                mu_rhk = np.poly1d(rhk_fits[i])(bv)
                sig_rhk = np.poly1d(scatter_fits[i])(bv)
                g = gaussian(mu_rhk,sig_rhk,RHK)
                plt.plot(RHK,g,label = CLUSTER_NAMES[i])


def plot_vs_age(y,title,y_label,axis,pp):
        plt.title(title)
        plt.semilogx(CLUSTER_AGES,y,'o')
        plt.xlabel(u'Age (Myr)',size=18)
        plt.ylabel(y_label,size=18)
        interp = np.interp(np.log10(CLUSTER_AGES),np.log10(CLUSTER_AGES),y)
        plt.semilogx(CLUSTER_AGES, interp, dashes=[2, 2])
        plt.axis(axis)
        #plt.semilogx(np.unique(ages), np.poly1d(np.polyfit(np.log10(ages), y, 1))(np.unique(np.log10(ages))), dashes=[2, 2])
        pp.savefig()
        #plt.close()
        plt.show()
"""

def plot_mamajek(bv_rhk,fits):
        import ca_constants as const
        #plt.figure(figsize=(7,6))
        plt.xlabel(r'$(B-V)_0$',size=18)
        #plt.ylabel(u'logR\'HK',size=18)
        #plt.ylabel(u'Log Calcium Abundance',size=17)
        plt.ylabel(r"Log(R'$_{HK})$",size=16)
        plt.axis([.45,.9,-5,-3.7])
        for i in range(len(bv_rhk)):
                color = const.COLORS[i]
                plt.scatter(bv_rhk[i][0],bv_rhk[i][1], color=color,marker=const.MARKERS[i],label=const.CLUSTER_NAMES[i])
                plt.plot(const.BV, fits[i][0](const.BV), color=color, dashes=[2, 2], alpha = .7)
                #shade_scatter(fits[i],const.BV.tolist())
        #plt.title('Calcium Abundance vs. B-V Color',size=18)
        plt.text(.7,-4.005,"Sco-Cen",size=13,color=const.COLORS[0])
        plt.text(.725,-4.33,"Pleiades",size=13,color=const.COLORS[1])
        plt.text(.75,-4.55,"Hyades",size=13,color=const.COLORS[2])
        plt.text(.61,-4.77,"M67",size=13,color=const.COLORS[3])


def fit_histogram(bv_m,fits,metal,pdfPage=None,showPlots=False,title=None,upper_limits=None,li_range=None,age_range=None):
    const = init_constants(metal)
    allClusters, totalStars = my_fits.get_fit_residuals(bv_m,fits,metal,upper_limits,li_range,age_range=age_range,scale_by_std= metal=='calcium')
    #allClusters, totalStars = my_fits.get_fit_residuals(bv_m,fits,metal,upper_limits,li_range,linSpace=True)
    
    pdf_fit,cdf_fit = my_fits.fit_histogram(metal,totalStars,fromFile=False,saveToFile=False)

    mu = np.mean(totalStars)
    sigma = np.std(totalStars)
    
    max_val = np.max([np.max(totalStars),np.abs(np.min(totalStars))]) + .1
    
    plt.rcParams["figure.figsize"] = (8,6)
    plt.hist(allClusters,bins=30,stacked=True,density=True,color=const.COLORS)
    #plt.hist(allClusters,bins=50,stacked=True,density=True)
   
    #plt.hist(np.power(10,bv_m[0][1]),bins=30,density=True)
    
    #for m in [1,1.5,2,2.5,3]:
    #    pdf = pdf_fit(np.log10(const.METAL)-m)[0,0,:]/const.METAL[0,0,:]
    #    prob.normalize(const.METAL[0,0,:],pdf)
    #    plt.plot(const.METAL[0,0,:],pdf,label='pdf')
    #    plt.plot(const.METAL[0,0,:],prob.normalize(const.METAL[0,0,:],prob.lognorm(const.METAL[0,0,:]/10**m,s=0.17)/10**m),label='lognorm')
        
        #g_pdf = prob.gaussian(np.log10(const.METAL[0,0,:])-m,0,0.17)/const.METAL[0,0,:]
        #prob.normalize(const.METAL[0,0,:],g_pdf)
        #plt.plot(const.METAL[0,0,:],g_pdf,label='shifted gauss')
    #    plt.xlim(0,min(3*10**m,1600))
    #    plt.legend()
    #    plt.show()

    x = np.linspace(-1*max_val,max_val,1000)
    #plt.plot(x,cdf_fit(x))
    plt.plot(x,pdf_fit(x),color='red',linewidth=2)
    #m = -3
    #clust_scatter = [fits[i][1](0.65) for i in range(len(fits))]
    #fit_gaussian = my_fits.fit_gaussian(np.log10(const.CLUSTER_AGES),clust_scatter)
    #plt.plot(x,pdf_fit((x - m)/fit_gaussian(10))/fit_gaussian(10),label='10')
    
    #plt.plot(x,pdf_fit((x-m)/fit_gaussian(100))/fit_gaussian(100),label='100')
    #plt.plot(x,pdf_fit((x-m)/fit_gaussian(1000))/fit_gaussian(1000),label='1000')
    #print fit_gaussian(10),fit_gaussian(100),fit_gaussian(1000)
    
    
    plt.plot(x,prob.gaussian(x,mu,sigma),linestyle = '--',color='gray')
    #plt.legend()
    #plt.show()

    #plt.legend([r'Gaussian ($\mu,\sigma$) = (%.2f,%.2f)' % (mu,sigma)] + const.CLUSTER_NAMES)
    #plt.legend([r'Gaussian ($\mu$=%.3f,$\sigma$=%.3f)' % (mu,sigma)] + const.CLUSTER_NAMES)
    plt.legend(['Astrophysical PDF'] + ['Best-fit Gaussian'] + const.CLUSTER_NAMES)
    plt.xlim([-max_val,max_val])
    x_axis = 'Log(EWLi) - Li Fit'
    if (metal.lower()[0] == 'c'):
        x_axis = "Log(R'" + r'$_{HK})$' + ' - Ca Fit'
    plt.xlabel(x_axis,size=18)
    plt.ylabel('Normalized Number of Stars',size=16)
    if (title):
        plt.title(title,size=18)
    plt.tight_layout()
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    #plt.close()        


def baffles_vs_mamajek(bv_rhk,fits,i,pdfPage=None,showPlots=False,title=None,mamaProduct=False):
    import ca_constants as const
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,fits,omit_cluster=i)
    my_ages = []
    my_error = []
    mamajek_ages = []
    post_prod = 0
    for j in range(len(bv_rhk[i][0])):
        b,r = bv_rhk[i][0][j], bv_rhk[i][1][j]
        mamajek_ages.append(utils.getMamaAge(r))
        post = baf.get_posterior(b,r)
        post_prod += np.log(post.array)
        stat = post.stats
        my_ages.append(stat[2])
        my_error.append((stat[2] - stat[1],stat[3] - stat[2]))
    
    post_prod = prob.normalize(const.AGE,np.exp(post_prod))
    baffles_age = prob.stats(const.AGE,post_prod)[2]
    
    plt.Line2D([0], [0], color='C%d'% i,marker=const.MARKERS[i],label=const.CLUSTER_NAMES[i])
    plt.axis([.4,14000,.4,14000])
    plt.title(const.CLUSTER_NAMES[i],size=16)
    plt.xlabel(r'Age derived from M & H (2008)',size=16)
    plt.ylabel(u'BAFFLES Age (Myr)',size=18)
    for j in range(len(my_ages)):
        if (j==0):
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,color=const.COLORS[i],marker=const.MARKERS[i])
        else:
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,color=const.COLORS[i],marker=const.MARKERS[i])
    
    plt.hlines(y=baffles_age,xmin= 0,xmax= baffles_age, \
        label='BAFFLES cluster age: %.3g Myr' % baffles_age ,linestyle='--',color = 'C0')
    

    if mamaProduct:
        age = utils.getMamaProductAge(bv_rhk[i][1])
        plt.vlines(x=age,ymin= 0,ymax= age, \
            label='M & H (2008): %.3g Myr' % age ,linestyle='--',color = 'C2')
    
    
    plt.plot([0,10000],[0,10000],dashes=[2,2],color='k')
    plt.plot(const.CLUSTER_AGES[i],const.CLUSTER_AGES[i],marker='*',markersize=18,color='C1',linestyle='None',label='Isochronal cluster age: %d Myr' % const.CLUSTER_AGES[i],alpha=1,zorder=10)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    #plt.legend(handles=legend_elements,loc=2)
    plt.legend(loc=2)
    plt.tight_layout()
    
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    plt.close()        
