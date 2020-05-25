"""
Adam Stanford-Moore
5/22/20
Functions creating useful plots for visualizing posteriors and fits
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy import interpolate
import bisect
import ca_constants as const
import probability as prob
import fitting as my_fits
import baffles
import utils
TITLE_SIZE = 16
AXIS_LABEL_SIZE = 16


def posterior(age,y,stat,title=' ',pp=None,showPlot=False,starArray = [],\
        givenAge=None,givenErr = None,mamajekAge=None,logPlot=False,bv_arr = None,\
        metal='calcium'):
    const = init_constants(metal)
    isUpperLim = len(stat) == 4 #arbitrarily age at CDF=[.002,.05,.34,1]
    plt.title(title,size=TITLE_SIZE)
    
    cmap,norm,sc = None,None,None
    if bv_arr is not None:
        cmap = plt.cm.get_cmap('RdYlBu_r')
        norm = mpl.colors.Normalize(vmin=const.BV_RANGE[0], vmax=const.BV_RANGE[1])
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
        print('Isochronal age exists within %f %% CI' % prob.get_percentile(age,y,givenAge))
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
    plt.xlabel(u'Age (Myr)',size=AXIS_LABEL_SIZE)
    plt.ylabel('Probability Density (Myr^-1)',size=AXIS_LABEL_SIZE)
    plt.legend() #loc="upper right")
    if (not logPlot and not isUpperLim):
        r = getAgeRange(stat,starArray,givenAge)
        plt.xlim(r)
    if (len(starArray) > 1):
            plt.ylim([0,np.max(y)*1.5])
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
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

#Determines the domain for plotting posterior in linear space
def getAgeRange(stat,starArray,givenAge):
        l = 0
        u = const.GALAXY_AGE
        median = stat[2]
        sigma = (stat[3] - stat[1])/2
        if (len(starArray) > 1):
            if (sigma < 5):
                sigma = 5
        if (median - 4*sigma > l):
            l = median - 4*sigma
        if (median + 4*sigma < u):
            u = median + 4*sigma
        if givenAge and not (l < givenAge < u):
            if givenAge < l:
                l = givenAge - 0.1*(median - l) #move lower-bound a little left
            elif givenAge > u:
                    u = givenAge + 0.1*(u - median) #move lower-bound a little left
        return [l,u]


def shade_scatter(fit,BV):
    y_top,y_bottom = fit[0](BV), fit[0](BV)
    y_top += fit[1](BV)
    y_bottom -= fit[1](BV)
    plt.fill_between(BV,y_bottom,y_top,alpha=.3)

def plot_fits(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,
              shadeScatter=True,specific_clusters=None,residuals=False):
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

        plt.plot(const.BV,fits[i][0](const.BV),dashes=[2,2], color=color)
        
        if shadeScatter:
            shade_scatter(fits[i],const.BV.tolist())
        plt.title(const.CLUSTER_NAMES[i],size=TITLE_SIZE)
        plt.axis(const.BV_RANGE + const.METAL_RANGE)
        plt.xlabel(r'$(B-V)_o$',size=AXIS_LABEL_SIZE)
        set_ylabel(metal)
        #plt.legend()
        plt.tight_layout()
        plt.minorticks_on()
        plt.tick_params(axis='both',which='both',right=True,top=True)
        if (pdfPage):
            pdfPage.savefig()
        if (showPlots):
            plt.show()
        plt.close()
        if residuals:
            resid = residuals(bv_m[i][0],bv_m[i][1],fits[i][s][0])
            plt.scatter(bv_m[i][0],resid,label = fit_names[s])
            plt.axhline(y=0,linestyle='--',color='k')
            plt.title(CLUSTER_NAMES[i] + ' Residuals',size=TITLE_SIZE)
            plt.xlim(const.BV_RANGE)
            plt.xlabel(r'$(B-V)_o$',size=AXIS_LABEL_SIZE)
            set_ylabel(metal)
            plt.legend()
            if (pdfPage):
                pdfPage.savefig()
            if (showPlots):
                plt.show()
            plt.close()


def plot_fits_subplot(plt,i,bv_m,fits,metal,upper_lim=None):
    const = init_constants(metal)   

    for j in range(len(bv_m[i][0])):
        if (upper_lim is not None and upper_lim[i][j]):
            plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=const.COLORS[i],marker=const.DOWN_ARROW)
        else:
            plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=const.COLORS[i],marker=const.MARKERS[i])

    plt.plot(const.BV,fits[i][0](const.BV),dashes=[2,2], color=const.COLORS[i])
    plt.scatter([],[],color=const.COLORS[i],marker=const.MARKERS[i],label=const.CLUSTER_NAMES[i])
    plt.legend()
    #plt.set_title(const.CLUSTER_NAMES[i],size=TITLE_SIZE)
    plt.axis(const.BV_RANGE + const.METAL_RANGE)

    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)



def metal_vs_bv(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,
                shadeScatter=False,title=None,primordial_li = False,fits_only=False,
                specific_clusters = None,legend=True,textlabels=False):
    const = init_constants(metal)
    plt.xlabel(r'$(B-V)_0$',size=AXIS_LABEL_SIZE)
    set_ylabel(metal)    
    plt.axis([const.BV_RANGE[0],const.BV_RANGE[1]] + const.METAL_RANGE)
    
    #plot primordial lithium
    if (primordial_li and metal.lower()[0]=='l'):
        #ngc2264_fit = fits[const.CLUSTER_NAMES.index("NGC2264")][0]
        primordial_li_fit = my_fits.MIST_primordial_li()#ngc2264_fit,fromFile=False,saveToFile=True) 
        plt.plot(const.BV,primordial_li_fit(const.BV),color=const.PRIM_COLOR,label="Primordial Li EW")
    
    clusters_to_plot = specific_clusters if specific_clusters != None else range(len(bv_m))
    for c in clusters_to_plot:
        ax = plt.gca()
        if (not fits_only):
            for i in range(len(bv_m[c][0])):
                if (upper_lim and upper_lim[c][i]):
                    plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=const.COLORS[c],
                                marker=const.DOWN_ARROW)
                else:
                    plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=const.COLORS[c],
                                marker=const.MARKERS[c])
        if legend:
            plt.scatter([],[],color=const.COLORS[c],marker=const.MARKERS[c],
                        label=const.CLUSTER_NAMES[c])
        
        if fits_only:
            plt.plot(const.BV,fits[c][0](const.BV), color=const.COLORS[c])
        else:
            plt.plot(const.BV,fits[c][0](const.BV),dashes=[2,2], color=const.COLORS[c])
        
        if (shadeScatter):
            shade_scatter(fits[c],const.BV.tolist())

    if textlabels and metal == 'calcium':
        plt.text(.7,-4.005,"Sco-Cen",size=13,color=const.COLORS[0])
        plt.text(.725,-4.27,"Pleiades",size=13,color=const.COLORS[5])
        plt.text(.75,-4.55,"Hyades",size=13,color=const.COLORS[7])
        plt.text(.61,-4.77,"M67",size=13,color=const.COLORS[8])

    if (title):
        plt.title(title,size=TITLE_SIZE)
    
    if legend:
        leg = plt.legend(loc='lower right')
        for i in range(len(leg.get_texts())):
            if (primordial_li and metal == 'lithium'):
                if (i==0):
                    plt.setp(leg.get_texts()[0],color=const.PRIM_COLOR)
                else:
                    plt.setp(leg.get_texts()[i],color=const.COLORS[clusters_to_plot[i-1]])
            else:
                plt.setp(leg.get_texts()[i],color=const.COLORS[clusters_to_plot[i]],size=13)
        
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
        
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()

#RHK vs age polynomial
def metal_vs_age_subplot(plt,fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,\
        shadeScatter=False,errorbars=True,bv_m=None,upper_lim=None, \
        mamajek_poly=False,metal_val=None,logAge = True,plotStars=False,omit_cluster=None,legend=True):
    const = init_constants(metal)
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = my_fits.get_valid_metal(bv,fits,const,None,None,omit_cluster)

    #plt.figure(figsize=(7,5))
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

    mu,sig,lbl = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,metal,omit_cluster)

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
        plt.fill_between(const.AGE,mu(const.AGE) - sig(const.AGE),
                mu(const.AGE) + sig(const.AGE),alpha=.2,color='C0')

    if (mamajek_poly and metal[0].lower()=='c'):
        plt.plot(const.AGE,utils.getMamaRHK(const.AGE),linestyle='--',color='gray',\
                label='M & H 2008')

    #ax = plt.gca()
    if logAge: plt.set_xscale('log')
    plt.axis([1,const.GALAXY_AGE]+const.METAL_RANGE)
    if legend:
        plt.legend()
    #plt.set_xlabel('Age (Myr)',size=AXIS_LABEL_SIZE)
    #ylabel = set_ylabel(metal)
    if (title):
        plt.set_title(title,size=TITLE_SIZE)
    #plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    #if (pdfPage):
    #    pdfPage.savefig()
    #if (showPlots):
    #    plt.show()
    #if pdfPage:
    #    plt.close()

#RHK vs age polynomial
def metal_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,\
        shadeScatter=False,errorbars=True,bv_m=None,upper_lim=None, \
        mamajek_poly=False,metal_val=None,logAge = True,plotStars=False,omit_cluster=None):
    const = init_constants(metal)
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = my_fits.get_valid_metal(bv,fits,const,None,None,omit_cluster)

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
    
    mu,sig,lbl = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,metal,omit_cluster)

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
        plt.fill_between(const.AGE,mu(const.AGE) - sig(const.AGE),
                mu(const.AGE) + sig(const.AGE),alpha=.2,color='C0')
    
    if (mamajek_poly and metal[0].lower()=='c'):
        plt.plot(const.AGE,utils.getMamaRHK(const.AGE),linestyle='--',color='gray',\
                label='M & H 2008')
    
    ax = plt.gca()
    if logAge: ax.set_xscale('log')
    plt.axis([1,const.GALAXY_AGE]+const.METAL_RANGE)
    plt.legend()
    plt.xlabel('Age (Myr)',size=AXIS_LABEL_SIZE)
    ylabel = set_ylabel(metal)
    if (title):
        plt.title(title,size=TITLE_SIZE)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()

def scatter_vs_bv(fits,metal,pdfPage=None,showPlots=False,title=None):
    const = init_constants(metal)
    plt.xlabel(r'$(B-V)_0$',size=AXIS_LABEL_SIZE)
    #plt.ylabel(r'Li EW (m$\AA$)',size=AXIS_LABEL_SIZE)
    m = 'Log(LiEW/m$\AA$)'
    plt.ylabel('Scatter in ' + m,size=AXIS_LABEL_SIZE)
    for c in range(len(fits)):
        ax = plt.gca()
        color = const.COLORS[c]
        plt.plot(const.BV,fits[c][1](const.BV),dashes=[2,2],label=const.CLUSTER_NAMES[c], color=color)
    plt.legend()
    if (title):
        plt.title(title,size=TITLE_SIZE)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()


def scatter_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,bv_m=None,upper_lim=None,errorbars=True):
    const = init_constants(metal)
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = my_fits.get_valid_metal(bv,fits,const)
    
    N = np.array(const.NUM_STARS)
    for i in range(len(scatter)):
        if errorbars:            
            err = scatter[i]/np.sqrt(2*N[i] - 2)
            plt.errorbar(CLUSTER_AGES[i], scatter[i],yerr=err,color=const.COLORS[i],\
                marker = const.MARKERS[i],markersize=8,capsize=4,zorder=10)
            plt.scatter([], [], color=const.COLORS[i],marker=const.MARKERS[i],label=CLUSTER_NAMES[i])
        else:
            plt.scatter(CLUSTER_AGES[i], scatter[i],color=const.COLORS[i],\
                marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    
    mu,sig,lbl = my_fits.vs_age_fits(bv,CLUSTER_AGES,rhk,scatter,metal)
   
    #lbl = 'gaussian fit' if metal=='calcium' else 'meas + astro' #r'constant $\sigma$ = %.3f' % fit(.65) 
    #plt.plot(const.AGE,sig(const.AGE),'--',label=lbl,color='orange')

    ax = plt.gca()
    ax.set_xscale('log')
    plt.legend()
    plt.xlabel('Age (Myr)',size=AXIS_LABEL_SIZE)
    y_axis  = 'Log(LiEW/m$\AA$)'
    if (metal.lower()[0] == 'c'):
        y_axis = "Log(R'" + r'$_{HK})$'
    plt.ylabel('Std. Dev. of ' + y_axis ,size=AXIS_LABEL_SIZE)
    if (title):
        plt.title(title,size=TITLE_SIZE)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()


def set_ylabel(metal):
    if (metal[0].lower() == 'c'):
        plt.ylabel("Log(R'" + r'$_{HK})$',size=AXIS_LABEL_SIZE)
        return 'Calcium Emission Strength'
    elif (metal[0].lower() == 'l'):
        plt.ylabel(r'Log(LiEW/m$\AA$)',size=AXIS_LABEL_SIZE)
        return 'Lithium Abundance'

def init_constants(metal):
    if (metal[0].lower() == 'c'):
        import ca_constants as const
    elif (metal[0].lower() == 'l'):
        import li_constants as const
    else:
        raise RuntimeError("No metal specified. Please enter lithium or calcium")
    return const


def fit_histogram(bv_m,fits,metal,pdfPage=None,showPlots=False,title=None,
                  upper_limits=None,li_range=None,age_range=None,plot_cdf=False):
    const = init_constants(metal)
    

    allClusters, totalStars = my_fits.get_fit_residuals(bv_m,fits,metal,upper_limits,
                            li_range,age_range=age_range,scale_by_std= False,vs_age_fit=True,zero_center=True)
    pdf_fit,cdf_fit = my_fits.fit_histogram(metal,totalStars,fromFile=False,saveToFile=False)
    #pdf_fit,cdf_fit = my_fits.fit_student_t(metal,totalStars,fromFile=False,saveToFile=True)
    
    mu = np.mean(totalStars)
    sigma = np.std(totalStars)
    
    max_val = np.max([np.max(totalStars),np.abs(np.min(totalStars))]) + .1
    x = np.linspace(-1*max_val,max_val,1000)

    if not plot_cdf:
        plt.plot(x,pdf_fit(x),color='red',linewidth=2)
        plt.plot(x,prob.gaussian(x,mu,sigma),linestyle = '--',color='gray')
        plt.hist(allClusters,bins=30,stacked=True,density=True,color=const.COLORS)
    else:        
        plt.plot(x,cdf_fit(x),color='red',linewidth=2,zorder=10)
        plt.plot(x,prob.gaussian_cdf(x,mu,sigma),linestyle = '--',color='gray',linewidth=2)
        cdf = np.array([(totalStars < n).sum() for n in x],dtype='float')
        cdf /= cdf[-1]
        
        plt.plot(x,cdf,linewidth=3)
        for c in range(len(allClusters)):
            cdf = np.array([(allClusters[c] < n).sum() for n in x],dtype='float')
            cdf /= cdf[-1]
            plt.plot(x,cdf,color=const.COLORS[c],alpha=.5)

    if plot_cdf: plt.legend(['Numerical Fit'] + ['Best-fit Gaussian'] + ["Total CDF"] \
                            + const.CLUSTER_NAMES)
    elif metal == "calcium": 
        plt.legend([r'Empirical ${\cal H}$(r|0)'] + ['Best-fit Gaussian'] + const.CLUSTER_NAMES)
    elif metal == "lithium": 
        plt.legend([r'Empirical ${\cal K}$($\ell$|0)'] + ['Best-fit Gaussian'] + const.CLUSTER_NAMES)

    plt.xlim([-max_val,max_val])
    #x_axis = 'Log(Li EW) - Li Fit'
    x_axis = 'Log(Li EW) - i(t,b)'
    y_axis = r'dp/dl (log(m$\AA$)^-1)'
    if (metal.lower()[0] == 'c'):
        #x_axis = "Log(R'" + r'$_{HK})$' + ' - Ca Fit'
        x_axis = "Log(R'" + r'$_{HK})$' + ' - f(t)'
        y_axis = 'Probability density'
    plt.xlabel(x_axis,size=AXIS_LABEL_SIZE)
    plt.ylabel(y_axis,size=AXIS_LABEL_SIZE)

    if plot_cdf: plt.ylabel('Cumulative probability',size=AXIS_LABEL_SIZE)
    if (title):
        plt.title(title,size=TITLE_SIZE)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()


def baffles_vs_mamajek(bv_rhk,fits,i,pdfPage=None,showPlots=False,title=None,mamaProduct=False):
    import ca_constants as const
    baf = baffles.age_estimator('calcium')
    #baf.make_grids(bv_rhk,fits,omit_cluster=i)
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
    plt.title(const.CLUSTER_NAMES[i],size=TITLE_SIZE)
    plt.xlabel(r'Age derived from M & H (2008)',size=AXIS_LABEL_SIZE)
    plt.ylabel(u'BAFFLES Age (Myr)',size=AXIS_LABEL_SIZE)
    for j in range(len(my_ages)):
        if (j==0):
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,
                         color=const.COLORS[i],marker=const.MARKERS[i])
        else:
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,
                         color=const.COLORS[i],marker=const.MARKERS[i])
    
    plt.hlines(y=baffles_age,xmin= 0,xmax= baffles_age, \
        label='BAFFLES cluster age: %.3g Myr' % baffles_age ,linestyle='--',color = 'C0')
    
    if mamaProduct:
        age = utils.getMamaProductAge(bv_rhk[i][1])
        plt.vlines(x=age,ymin= 0,ymax= age, \
            label='M & H (2008): %.3g Myr' % age ,linestyle='--',color = 'C2')
        
    plt.plot([0,10000],[0,10000],dashes=[2,2],color='k')
    plt.plot(const.CLUSTER_AGES[i],const.CLUSTER_AGES[i],marker='*',markersize=18,
             color='C1',linestyle='None',label='Isochronal cluster age: %d Myr'
             % const.CLUSTER_AGES[i],alpha=1,zorder=10)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.legend(loc=2)
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params(axis='both',which='both',right=True,top=True)
    
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    if pdfPage:
        plt.close()

def plot_mamajek(bv_rhk,fits):
    import ca_constants as const
    #plt.figure(figsize=(7,6))
    plt.xlabel(r'$(B-V)_0$',size=AXIS_LABEL_SIZE)
    #plt.ylabel(u'logR\'HK',size=AXIS_LABEL_SIZE)
    #plt.ylabel(u'Log Calcium Abundance',size=AXIS_LABEL_SIZE)
    plt.ylabel(r"Log(R'$_{HK})$",size=AXIS_LABEL_SIZE)
    plt.axis([.45,.9,-5,-3.7])
    for i in range(len(bv_rhk)):
        color = const.COLORS[i]
        plt.scatter(bv_rhk[i][0],bv_rhk[i][1], color=color,marker=const.MARKERS[i],label=const.CLUSTER_NAMES[i])
        plt.plot(const.BV, fits[i][0](const.BV), color=color, dashes=[2, 2], alpha = .7)
    plt.text(.7,-4.005,"Sco-Cen",size=13,color=const.COLORS[0])
    plt.text(.725,-4.33,"Pleiades",size=13,color=const.COLORS[1])
    plt.text(.75,-4.55,"Hyades",size=13,color=const.COLORS[2])
    plt.text(.61,-4.77,"M67",size=13,color=const.COLORS[3])


