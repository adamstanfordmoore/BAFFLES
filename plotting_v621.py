import matplotlib.pyplot as plt
import bisect
import ca_constants as const
import numpy as np
import probability as prob
import fitting as my_fits
from scipy import interpolate
import baffles
#COLORS = ['C0','C4','C2','C1','C5','C6','C7','C8']


#givenErr = [-5,5] or 5
def posterior(age,y,stat,title='Posterior Plot',pp=None,showPlot=False,starArray = [],givenAge=None,givenErr = None,mamajekAge=None,logPlot=False):
    isUpperLim = len(stat) == 4 #arbitrarily [.34,.45,.002,1]
    plt.title(title)
    for post in starArray:
        prob.scale_to_height(post,np.max(y))
        plt.plot(const.AGE,post,alpha = .4,linewidth=1,color='.5')
    
    if (givenAge):
        plt.axvline(x=givenAge,color='r',label='True Age: %d Myr' % givenAge,zorder=10)
        if (givenErr):
            if (type(givenErr) == float or type(givenErr) == int):
                givenErr = [-1*givenErr,givenErr]
            plt.axvspan(givenAge+givenErr[0], givenAge+givenErr[1], alpha=0.2, color='r')
    if (mamajekAge):
        plt.axvline(x=mamajekAge,color='C2',label='MH08 age: %d' % mamajekAge,zorder=10)
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
    plt.close()

def shadeStats(age,y,stat,upperLim):
    if upperLim:
        plt.fill_between(age,y, where= (age >= stat[0]) & (age <= stat[-1]),\
                color='.3', label=r'1$\sigma$ lower lim: %.2d Myr' % stat[0])
        plt.fill_between(age,y, where= (age >= stat[1]) & (age <= stat[0]),\
                color='.6',alpha=1, label=r'2$\sigma$ lower lim: %.2d Myr' % stat[1])
        plt.fill_between(age,y, where= (age >= stat[2]) & (age <= stat[1]),\
                color='.9',alpha=1,label=r'3$\sigma$ lower lim: %.2d Myr' % stat[2])
        return

    plt.vlines(x=stat[2],ymin= 0,ymax= y[bisect.bisect_left(const.AGE,stat[2])], \
            label='BAFFLES Age: %.2d Myr' % stat[2] ,color = 'orange')
    plt.fill_between(age,y, where= (age >= stat[1]) & (age <= stat[3]),color='.3', \
            label='68%% CI: %.2d - %.2d' % (stat[1],stat[-2]))
    plt.fill_between(age,y, where= (age >= stat[0]) & (age <= stat[-1]),color='.6',\
            alpha=0.5)#, label='95%% CI: %.3g - %.3g' % (stat[0],stat[-1]),alpha = .5)

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

def plot_fits(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,titles=None):
    if (not pdfPage and not showPlots):
        return
    const = init_constants(metal)   
    for i in range(len(fits)):
        ax = plt.gca()
        color = next(ax._get_lines.prop_cycler)['color']
        for j in range(len(bv_m[i][0])):
            if (upper_lim and upper_lim[i][j]):
                plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=color,marker=const.DOWN_ARROW)
            else:
                plt.scatter(bv_m[i][0][j],bv_m[i][1][j],color=color)

        #l = 'piecewise: %d' % SEG if converged[i][s] else 'Not Converged: %d' % SEG
        #l = 'piecewise, bin size: %d' % SEG if converged[i][s] else 'Not Converged, bin size: %d' % SEG
        plt.plot(const.BV,fits[i][0](const.BV),dashes=[2,2], color=color)
        shade_scatter(fits[i],const.BV.tolist())
        if (titles):
            plt.title(titles[i])
        else:
            plt.title(const.CLUSTER_NAMES[i] + r' Fit with Shaded $1\sigma$ Std Dev')
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

def metal_vs_bv(bv_m,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,shadeScatter=False,title=None,primordial_li = False,fits_only=False,specific_clusters = None):
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
        ngc2264_fit = fits[const.CLUSTER_NAMES.index("NGC2264")][0]
        primordial_li_fit = my_fits.primordial_li(ngc2264_fit,fromFile=False,saveToFile=True) 
        plt.plot(const.BV,primordial_li_fit(const.BV),label="Primordial LiEW")
    
    clusters_to_plot = specific_clusters if specific_clusters != None else range(len(bv_m))
    for c in clusters_to_plot:
            ax = plt.gca()
            color = const.COLORS[c] #next(ax._get_lines.prop_cycler)['color']
            marker = const.MARKERS[c]
            if (not fits_only):
                for i in range(len(bv_m[c][0])):
                    if (upper_lim and upper_lim[c][i]):
                        plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=color,marker=const.DOWN_ARROW)
                    else:
                        plt.scatter(bv_m[c][0][i],bv_m[c][1][i],color=color,marker=marker)
            plt.plot(const.BV,fits[c][0](const.BV),dashes=[2,2],label=const.CLUSTER_NAMES[c], color=color)
            
            if (c == 0):
                plt.text(.7,-4.005,"Sco-Cen",size=13,color=const.COLORS[0])
            elif c == 4:
                plt.text(.725,-4.33,"Pleiades",size=13,color=const.COLORS[4])
            elif c == 6:
                plt.text(.75,-4.55,"Hyades",size=13,color=const.COLORS[6])
            elif c == 7:    
                plt.text(.61,-4.77,"M67",size=13,color=const.COLORS[7])
    if (shadeScatter):
        shade_scatter(fits[c],const.BV.tolist())

    if (title):
        plt.title(title,size=18)
        #else:
        #   plt.title(metal[0].upper() + metal[1:] + r' Fits')
    
    if metal.lower()[0] == 'l':
        leg = plt.legend(loc=7)
        for i in range(len(leg.get_texts())):
            if (primordial_li):
                if (i==0):
                    plt.setp(leg.get_texts()[0],color='C0')
                else:
                    plt.setp(leg.get_texts()[i],color=const.COLORS[i-1])
            else:
                color = const.COLORS[clusters_to_plot[i]]
                plt.setp(leg.get_texts()[i],color=color)
    plt.tight_layout()
        
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    #plt.close()

#RHK vs age polynomial
def metal_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,shadeScatter=False,errorbars=False,bv_m=None,upper_lim=None, mamajek_poly=False):
    const = init_constants(metal)
    #rhk = [fits[i][0](bv) for i in range(len(fits))]
    #sig = np.mean([fits[i][1](bv) for i in range(len(fits))]) #Change if non-constant scatter
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = [],[],[],[]
    
    if (metal[0].lower()=='l'):
        #Add primordial lihtium point
        scatter.append(0)
        rhk.append(my_fits.primordial_li()(bv))
        CLUSTER_AGES.append(const.PRIMORDIAL_LI_AGE)
        CLUSTER_NAMES.append("Primordial LiEW")

    for i in range(len(fits)):
        r = fits[i][0](bv)
        if (const.METAL_RANGE[0] <  r < const.METAL_RANGE[1]):
            rhk.append(r)
            scatter.append(fits[i][1](bv))
            CLUSTER_AGES.append(const.CLUSTER_AGES[i])
            CLUSTER_NAMES.append(const.CLUSTER_NAMES[i])
    #sig = np.mean(scatter)
    if (metal[0].lower()=='l'):#info added from depletion boundary
        #add bldb point
        scatter.append(0) 
        rhk.append(const.ZERO_LI)
        CLUSTER_AGES.append(my_fits.bldb_fit(fits)(bv))
        CLUSTER_NAMES.append('BLDB point')

    if (metal[0].lower()=='l'):
        plt.figure(figsize=(7,5))
    ax = plt.gca()
    ax.set_xscale('log')
    plt.ylim(const.METAL_RANGE)
    for i in range(len(rhk)):
        if (errorbars):
            plt.errorbar(CLUSTER_AGES[i], rhk[i], yerr=scatter[i],color=const.COLORS[i],marker='o',label=CLUSTER_NAMES[i],zorder=10)
        else:
            plt.scatter(CLUSTER_AGES[i], rhk[i],color=const.COLORS[i],marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    #for c in range(len(rhk)):
        #plt.errorbar(CLUSTER_AGES[c],rhk[c],yerr=scatter[c],marker=const.MARKERS[c],color=const.COLORS[c],label=CLUSTER_NAMES[c])
    
    fit,lbl,sig = None,None,None
    if (metal[0].lower() == 'c'):
        fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,2)
        lbl='Polynomial fit'
        sig = my_fits.total_scatter(bv_m,fits)
    elif (metal[0].lower() == 'l'):
        if (.76 < bv < .94):
            fit = my_fits.piecewise([0,2.2,3],[2.5,2.1,.5])
            #fit = my_fits.pwise_fit(np.log10(CLUSTER_AGES),rhk,segments=2,guess_fit=guess,x_method='free')[0]
            lbl = 'Piecewise'
        else:
            fit = my_fits.constrained_poly_fit(np.log10(CLUSTER_AGES),rhk,0)
            lbl = 'Constrained Polynomial'
        sig = my_fits.fit_two_scatters(bv_m,fits,upper_lim)(fit(np.log10(const.AGE)))
    
    
    plt.plot(const.AGE,fit(np.log10(const.AGE)),label=lbl)
    """ 
    if (len(rhk) < 3):
        fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,len(rhk) - 1)
        plt.plot(const.AGE,fit(np.log10(const.AGE)))
    elif (metal[0].lower() == 'c'):
        fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,2)
        plt.plot(const.AGE,fit(np.log10(const.AGE)),label='Polynomial fit')
    elif (metal[0].lower() == 'l'):
        #fit = interpolate.interp1d(np.log10(CLUSTER_AGES),rhk, fill_value='extrapolate')
        #plt.plot(const.AGE,fit(np.log10(const.AGE)),label='Interpolation')
        lbl = "linear interpolate"
        if (bv <=1.35):
            fit = interpolate.interp1d(np.log10(CLUSTER_AGES),rhk, fill_value='extrapolate')
        elif (bv <= 1.55):    
            fit = my_fits.constrained_poly_fit(np.log10(CLUSTER_AGES),rhk,0)
            lbl = 'constrained poly fit' 
        else:
            fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),rhk,1)
            lbl = 'linear fit' 
        plt.plot(const.AGE,fit(np.log10(const.AGE)),label=lbl)
    """
    if (shadeScatter):
        plt.fill_between(const.AGE,fit(np.log10(const.AGE)) - sig,fit(np.log10(const.AGE)) + sig,alpha=.2,color='C0')
        #plt.scatter(CLUSTER_AGES[i], rhk[i],color=const.COLORS[i],marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    
    if (mamajek_poly and metal[0].lower()=='c'):
        plt.plot(const.AGE,my_fits.getMamaRHK(const.AGE),linestyle='--',color='gray',label='Mamajek & Hillenbrand 2010')
    
    
    
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
    m = 'Log(Li EW/m$\AA$)'
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
    #rhk = [fits[i][0](bv) for i in range(len(fits))]
    #sig = [fits[i][1](bv) for i in range(len(fits))])] 
    ax = plt.gca()
    ax.set_xscale('log')
    rhk,sig,CLUSTER_AGES,CLUSTER_NAMES = [],[],[],[]
    for i in range(len(fits)):
        r = fits[i][0](bv)
        if (const.METAL_RANGE[0] <  r < const.METAL_RANGE[1]):
            rhk.append(r)
            sig.append(fits[i][1](bv))
            CLUSTER_AGES.append(const.CLUSTER_AGES[i])
            CLUSTER_NAMES.append(const.CLUSTER_NAMES[i])
            print(rhk[-1],CLUSTER_NAMES[-1])

    for i in range(len(sig)):
        plt.scatter(CLUSTER_AGES[i], sig[i],color=const.COLORS[i],marker = const.MARKERS[i],label=CLUSTER_NAMES[i],s=80,zorder=10)
    fit,lbl = None,None
    if metal == 'calcium':
        fit = my_fits.fit_gaussian(np.log10(CLUSTER_AGES),sig)
        lbl = 'gaussian fit'
        plt.plot(const.AGE,fit(np.log10(const.AGE)),'--',label=lbl,color='orange')
    elif metal=='lithium':
        assert type(bv_m) != type(None) and type(upper_lim) != type(None), "Need bv_m and upper_lim"
        #fit = my_fits.poly_fit(np.log10(CLUSTER_AGES),sig,2)
        # remember function of EW not B-V
        print bv
        EW_fit = None
        if (.76 < bv < .94):
            EW_fit = my_fits.piecewise([0,2.2,3],[2.5,2.1,.5])
        else:
            EW_fit = my_fits.constrained_poly_fit(np.log10(CLUSTER_AGES),rhk,0)
    
        plt.plot(const.AGE,EW_fit(np.log10(const.AGE)))
        
        fit = my_fits.fit_two_scatters(bv_m,fits,upper_lim)(EW_fit(np.log10(const.AGE)))
        lbl = 'meas + astro' #r'constant $\sigma$ = %.3f' % fit(.65) 
        #plt.plot(const.AGE,fit,'--',label=lbl,color='orange')


    plt.legend()
    plt.xlabel('Age (Myr)',size=18)
    x_axis  = 'Log(Li EW/m$\AA$)'
    if (metal.lower()[0] == 'c'):
        x_axis = 'Log(R\'HK)'
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
        plt.ylabel(u'Log(R\'HK) -- Ca Emission Strength',size=16)
        return 'Calcium Emission Strength'
    elif (metal[0].lower() == 'l'):
        #plt.ylabel('Lithium Abundance',size=16)
        #plt.ylabel(u'log(EW/mA)',size=18)
        plt.ylabel(r'Log(Li EW/m$\AA$)',size=18)
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


def fit_histogram(bv_m,fits,metal,pdfPage=None,showPlots=False,title=None,upper_limits=None,li_range=None):
    const = init_constants(metal)
    allClusters = []
    #totalStars = []

    for c in range(len(fits)):
        arr = []
        residuals = my_fits.residuals(bv_m[c][0],bv_m[c][1],fits[c][0])
        for i in range(len(residuals)):
            if (upper_limits and upper_limits[c][i]):
                continue
            if (li_range and (bv_m[c][1][i] < li_range[0] or li_range[1] < bv_m[c][1][i])):
                continue
            arr.append(residuals[i])
        allClusters.append(arr)

    totalStars = np.concatenate(allClusters)
    mu = np.mean(totalStars)
    sigma = np.std(totalStars)
    
    max_val = np.max([np.max(totalStars),np.abs(np.min(totalStars))]) + .1
    
    plt.rcParams["figure.figsize"] = (8,6)
    plt.hist(allClusters,bins=20,stacked=True,density=True)
    
    x = np.linspace(-1*max_val,max_val,100)
    plt.plot(x,prob.gaussian(x,mu,sigma))

    #plt.legend([r'Gaussian ($\mu,\sigma$) = (%.2f,%.2f)' % (mu,sigma)] + const.CLUSTER_NAMES)
    plt.legend([r'Gaussian ($\mu$=%.3f,$\sigma$=%.3f)' % (mu,sigma)] + const.CLUSTER_NAMES)
    plt.xlim([-max_val,max_val])
    x_axis = 'Log(Li EW) - Li Fit'
    if (metal.lower()[0] == 'c'):
        x_axis = 'Log(R\'HK) - Ca Fit'
    plt.xlabel(x_axis,size=18)
    plt.ylabel('Number of Stars, Normalized',size=16)
    if (title):
        plt.title(title,size=18)
    plt.tight_layout()
    if (pdfPage):
        pdfPage.savefig()
    if (showPlots):
        plt.show()
    plt.close()        


def baffles_vs_mamajek(bv_rhk,fits,i,pdfPage=None,showPlots=False,title=None):
    import ca_constants as const
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_rhk,fits,omit_cluster=i)
    my_ages = []
    my_error = []
    mamajek_ages = []
    for j in range(len(bv_rhk[i][0])):
        b,r = bv_rhk[i][0][j], bv_rhk[i][1][j]
        mamajek_ages.append(my_fits.getMamaAge(r))
        post = baf.get_posterior(b,r)
        stat = post.stats
        my_ages.append(stat[2])
        my_error.append((stat[2] - stat[1],stat[3] - stat[2]))

    plt.Line2D([0], [0], color='C%d'% i,marker=const.MARKERS[i],label=const.CLUSTER_NAMES[i])
    plt.axis([.4,14000,.4,14000])
    plt.title(const.CLUSTER_NAMES[i] + " Stellar Age Comparison",size=18)
    plt.xlabel(r'Non-Bayesian Age (Myr)',size=18)
    plt.ylabel(u'BAFFLES Age (Myr)',size=18)
    for j in range(len(my_ages)):
        if (j==0):
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,label=const.CLUSTER_NAMES[i] + ' stars',fmt=const.FMTS[i])
        else:
            plt.errorbar(mamajek_ages[j],my_ages[j],np.array([my_error[j]]).T,fmt=const.FMTS[i])
    plt.plot([0,10000],[0,10000],dashes=[2,2],color='k')
    plt.plot(const.CLUSTER_AGES[i],const.CLUSTER_AGES[i],marker='*',markersize=18,color='C1',linestyle='None',label='Accepted Age: %d' % const.CLUSTER_AGES[i],alpha=1,zorder=10)
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
