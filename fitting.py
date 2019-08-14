"""
Module giving necessary fitting functions
"""

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate,integrate
from scipy.stats import norm
from scipy.optimize import minimize,curve_fit
import copy
import bisect
import pickle
import probability as prob
from astropy.io import ascii
from astropy.table import Table
import utils
import readData
from scipy.signal import savgol_filter
import time

BIN_SIZE = 10
MIN_PER_BIN = 4
MIN_LENGTH = .15 #minimum B-V length of piecewise segment
DOWN_ARROW = u'$\u2193$'
MEASUREMENT_LI_SCATTER = 10 #mA in linear space from measurement error
PRIMORDIAL_NLI = 3.2

#limits vertex from moving to the right of lim
# or sets vertex to be x1,y1. sets vertex to be first point by default
def constrained_poly_fit(x,y,x1=None,y1=None,lim=None,sigma=None):
    if x1 is None and lim is None: x1 = x[0]
    if y1 is None and lim is None: y1 = y[0]
    guess = np.polyfit(x,y,2)[0] if lim is None else np.polyfit(x,y,2)
    res = minimize(constrained_poly_minimizer,guess, args=(x,y,x1,y1,lim,sigma),method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing of constrained polynomial function, check initial guess"
        print res.x
    #print res.x
    return poly_constr_vert(res.x,x1,y1) if lim is None else np.poly1d(res.x)

def right_cubic_root(params):
    root = 4*params[1]*params[1] - 12*params[0]*params[2]
    if (root < 0):
        return -100000
        #raise RuntimeError("CUBIC ERROR")
    a = (-2*params[1] - np.sqrt(root)) / (6 * params[0])
    b = (-2*params[1] + np.sqrt(root)) / (6 * params[0])
    if (a >= b):
        return a
    return b
    
def poly_constr_vert(param,x1,y1):
    a = param[0]
    b = -2*a*x1
    c = y1 - a*x1**2 - b*x1
    return np.poly1d([a,b,c])

def constrained_poly_minimizer(params,x,y,x1,y1,lim,sigma):
    y_model = None
    if lim is not None:
        #if max is greater than lim return a large number
        vert = -params[1]/(2*params[0])
        if (vert > lim):
            return 1e10*(vert - lim)
        y_model = np.poly1d(params)(x)
    else:
        y_model = poly_constr_vert(params,x1,y1)(x) #np.poly1d(params)(x)
    return prob.chi_sqr(y,y_model,sigma,total=True)

"""
def constrained_poly_minimizer(params,x,y,scatter,lim):
    #if max is greater than lim return a large number
    vert = -params[1]/(2*params[0])
    #vert = right_cubic_root(params)
    if (vert > lim):
        return 1e10*(vert - lim)
    y_model = np.poly1d(params)(x)
    return prob.chi_sqr(y,y_model,scatter,total=True)
"""

#computes fit of age and lithium depletion boundary (well bv at which lithium goes to zero)
# returns function such that giving it a b-v value returns oldest age it could be.  
def bldb_fit(fits,plot=False): 
    bv_at_zero_li,ages,cluster_names = [],[],[]
    import li_constants as const
    
    BV = np.linspace(0.6,2.2,300) #picked to include M35
    for c in range(1,len(fits)): #omitting NGC2264
        li = fits[c][0](BV)
        i = np.argmin(np.abs(li - const.ZERO_LI))
        if (i != len(li)):
            bv_at_zero_li.append(BV[i])
            ages.append(const.CLUSTER_AGES[c])
            cluster_names.append(const.CLUSTER_NAMES[c])
   
    #print cluster_names
    upper_clusters = ['M67','Hyades','M34','M35']
    upper_ages = [ages[cluster_names.index(name)] for name in upper_clusters]
    upper_bv = [bv_at_zero_li[cluster_names.index(name)] for name in upper_clusters]
    fit = piecewise(upper_bv,np.log10(upper_ages))
    
    #fit = poly_fit(bv_at_zero_li,np.log10(ages),1)
    #fit = pwise_fit(bv_at_zero_li,np.log10(ages),2,guess_fit=linear_fit,x_method='free')[0]
    
    ldb_age = lambda x : np.power(10,fit(x))
    
    if (plot):
        pp = PdfPages('bldb_vs_age.pdf')
        plt.xlabel('B-V',size=18)
        plt.ylabel('Age (Myr)',size=18)
        plt.yscale('log')
        for c in range(len(ages)):
            plt.scatter(bv_at_zero_li[c],ages[c],s=60,label=cluster_names[c],color=const.COLORS[c+1],marker=const.MARKERS[c+1])
        plt.plot(BV,ldb_age(BV))
        plt.fill_between(BV,ldb_age(BV),color='C0',alpha=.2,label="valid ages")
        plt.legend()
        pp.savefig()
        plt.show()
        pp.close()
    
    return ldb_age

# not sure this works, finds scatter in age instead of Li ?
#def ldb_scatter(fits):
#    bv_at_zero_li,ages,cluster_names = [],[],[]
#    import li_constants as const
#    for c in range(len(fits)):
#        li = fits[c][0](const.BV)
#        i = bisect.bisect_left(const.BV,.65) #only interested in zero crossing at bv > .65
#        while (i < len(li) and li[i] > const.ZERO_LI):
#            i += 1 #can change to 5 for slightly more speed
#        if (i != len(li)):
#            bv_at_zero_li.append(const.BV[i])
#            ages.append(const.CLUSTER_AGES[c])
#            cluster_names.append(const.CLUSTER_NAMES[c])
#    fit = poly_fit(bv_at_zero_li,np.log10(ages),1)
#    def ldb_age(x):
#        return np.power(10,fit(x))
#    return np.std(residuals(bv_at_zero_li,np.log10(ages),ldb_age))

#returns the 1d polynomial
def linear_fit(x,y,scatter=False):
    return poly_fit(x,y,1,scatter=scatter)

def dip_gaussian(x,mu,sig,A,C):
    chi2 = prob.chi_sqr(x,mu,sig)
    return C + A/(sig*np.sqrt(2*np.pi)) * np.exp(-chi2/2)

def dip_poly(x,a,b,c):
    return a*x**2 + b*x + c

#def li_dip_fit(bv,li,upper_lim,dip_bv_range,dip_li_range):
#    const = utils.init_constants('lithium')
#    # separate out dip stars from non dip
#
#    bv,li,upper_lim = np.array(bv),np.array(li),np.array(upper_lim)
#    dip_mask = (bv <= dip_bv_range[1]) & (bv >= dip_bv_range[0]) & \
#               (li <= dip_li_range[1]) & (li >= dip_li_range[0])
#    mask = np.invert(dip_mask)
#    bv_dip,li_dip,upper_lim_dip = bv[dip_mask],li[dip_mask],upper_lim[dip_mask]
#    bv_,li_,upper_lim_ = bv[mask],li[mask],upper_lim[mask]
#    no_dip_fit,no_dip_sig_fit = poly_fit(bv_,li_,2,upper_lim_)
#    resid = residuals(bv,li,no_dip_fit)
#    popt,pcov = curve_fit(dip_gaussian,bv,resid,p0=[0.45,.04,-.1,0])
#    plt.scatter(bv,resid)
#    xx = np.linspace(min(bv),max(bv),100)
#    plt.plot(xx,dip_gaussian(xx,*popt))
#    plt.show()
#
#    def dip_fit(x): 
#        return no_dip_fit(x) + dip_gaussian(np.array(x),*popt)
#    sig = np.std(residuals(bv,li,dip_fit))
#    plt.scatter(bv,li)
#    plt.plot(const.BV,dip_fit(const.BV))
#    plt.show()

# dip_bv_range and dip_li_range specify box containing dipped stars,
# so they can be excluded for background polynomial
# total_dip_fit_range specifies bv range to fit polynomial
## edge_box specifies box of stars (bv_min,bv_max,li_min,li_max) to make dip go through
def li_dip_fit(bv,li,upper_lim,dip_bv_range,dip_li_range,dip_fit_range,edge_box):
    const = utils.init_constants('lithium')
    # separate out dip stars from non dip

    bv,li,upper_lim = np.array(bv),np.array(li),np.array(upper_lim)
    dip_mask = (bv <= dip_bv_range[1]) & (bv >= dip_bv_range[0]) & \
               (li <= dip_li_range[1]) & (li >= dip_li_range[0])
    mask = np.invert(dip_mask)
    dip_mask = (dip_mask) | ((bv <= dip_fit_range[1]) & (bv >= dip_fit_range[0]))
    
    bv_dip,li_dip,upper_lim_dip = bv[dip_mask],li[dip_mask],upper_lim[dip_mask]
    bv_,li_,upper_lim_ = bv[mask],li[mask],upper_lim[mask]

    #plt.scatter(bv_dip,li_dip)
    #plt.show()
    #plt.scatter(bv_,li_)
    #plt.show()


    no_dip_fit,no_dip_sig_fit = poly_fit(bv_,li_,2,upper_lim_)
    
    edge_mask = (bv_dip > edge_box[0]) & (bv_dip < edge_box[1]) & \
            (li_dip > edge_box[2]) & (li_dip < edge_box[3])
    #edge_mask = (bv_dip < .507) & (bv_dip > 0.39) & (li_dip > 1.93) & (li_dip < 2.075) 

    sigma = np.where(edge_mask,.02,1) # makes fit go through edge_box points

    #dip_fit,dip_sig_fit = linear_fit(bv_dip,li_dip,scatter=True)
    popt,pcov = curve_fit(dip_poly,bv_dip,li_dip,sigma=sigma)#,p0=[100,-90,21.25])
    #print popt
    #popt,pcov = curve_fit(dip_gaussian,bv_dip,li_dip,p0=[0.45,.04,-.1],sigma=sigma)
    dip_fit = lambda x: dip_poly(np.array(x),*popt)
    #dip_fit = lambda x: dip_gaussian(np.array(x),*popt)
    dip_sig_fit = np.poly1d([.25])

    final_fit = lambda x : np.where((np.array(x) <= dip_bv_range[1]) & \
            (np.array(x) >= dip_bv_range[0]),dip_fit(x),no_dip_fit(x))
    
    final_sig_fit = lambda x: np.where((np.array(x) <= dip_bv_range[1]) & \
            (np.array(x) >= dip_bv_range[0]),dip_sig_fit(x),no_dip_sig_fit(x))
    
    return [final_fit,final_sig_fit]


# finds and returns n dimensional polynomial fit.  If scatter=True than it returns [median fit,scatter fit]
def poly_fit(x,y,n=2,upper_lim=None,scatter=False,weight = None):
    if (upper_lim is not None):
        return minimize_polynomial(x,y,n,upper_lim)
    fit = np.poly1d(np.polyfit(x,y,n,w=weight))
    if (scatter):
        scatter_fit = np.poly1d(np.std(residuals(x,y,fit)))
        return [fit,scatter_fit]
    return fit

def minimize_polynomial(x,y,n,upper_lim):
    guess_poly = poly_fit(x,y,n,scatter=True)
    guess = np.append(guess_poly[0],guess_poly[1])  
    res = minimize(poly_minimizer,guess,args=(x,y,upper_lim), method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing of polynomial function, check initial guess"
        print res.x

    #comp.append(compare_fits(res.x,argv))
    #plot_comp(comp,argv[3],res.success)

    sig = res.x[-1]
    poly = res.x[0:-1]
    fit = np.poly1d(poly)
    return [fit, constant_fit(sig)[0]]
    #return [fit, log_scatter_with_linear_offset(sig,np.poly1d(res.x[0:-1]))]


def poly_minimizer(params,c,r,upper_lim):
    sig = params[-1]
    if (utils.negative_sig(sig)):
        return np.inf
    r_model = np.poly1d(params[0:-1])(c)
    #sig_model = log_scatter_with_linear_offset(sig,np.poly1d(params[0:-1]))(c)
    sig_model = constant_fit(sig)[0](c)
       
    return inverse_log_likelihood(r_model,sig_model,r,upper_lim)

# 4 params = mu,sig,amplitude,vertical_offset
# takes in log ages and returns function of normal ages
def fit_gaussian(x,y):
    guess = [np.mean(x),np.std(x),np.max(y)-np.min(y),np.min(y)]
    res = minimize(gaussian_scatter_minimizer,guess,args=(x,y), method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing of gaussian scatter_vs_age fit"
        print res.x

    [mu,sig,A,c] = res.x
    print res.x
    def gauss_fit(x_coords):
        return prob.gaussian(np.log10(x_coords),mu,sig)*A + c
    
    if not res.success:
        import ca_constants as const
        plt.semilogx(const.AGE,gauss_fit(np.log10(const.AGE)))
        plt.show()
    
    return gauss_fit

def gaussian_scatter_minimizer(params,x,y):
    import ca_constants as const
    [mu,sig,A,c] = params
    fit = prob.gaussian(x,mu,sig)*A + c
    num_stars = np.array(const.CLUSTER_INDEX)[:,1] - np.array(const.CLUSTER_INDEX)[:,0]
    #print num_stars
    sigma_weight = 1.0/num_stars
    #print sigma_weight
    return prob.chi_sqr(y,fit,sigma_weight,total=True)

def constant_fit(y):
    m,s = np.mean(y),np.std(y)
    return [np.poly1d([m]),np.poly1d([s])]



# finds total scatter in log space combining astrophysical scatter and constant measurement uncertainty
# sig is the constant astrophysical scatter in log space, li_fit is the log lithium mean functional fit
# returns a function of log total scatter as a function of B-V
def log_scatter_with_linear_offset(sig,li_fit):
    import li_constants as const
    logEW = li_fit(const.BV)
    linEW = np.power(10,li_fit(const.BV))
    sig_a_lin = (np.power(10,logEW + sig) - np.power(10,logEW - sig))/2
    sig_tot_lin = np.sqrt(sig_a_lin**2 + MEASUREMENT_LI_SCATTER**2) 
    #upper_log = np.log10(linEW + sig_tot_lin)
    #lower = linEW - sig_tot_lin
    #lower[lower <=0] = 1
    sig_tot_log = np.log10(linEW + sig_tot_lin)
    return interpolate.interp1d(const.BV,sig_tot_log, fill_value='extrapolate')

#returns a single number for the average scatter across all the clusters
def total_scatter(bv_m,fits,omit_cluster=None,upper_limits=None,li_range=None,scale_by_scatter=None):
    allClusters = []
    for c in range(len(fits)):
        if (omit_cluster and c==omit_cluster):        
            continue    
        arr = []
        resid = residuals(bv_m[c][0],bv_m[c][1],fits[c][0])
        for i in range(len(resid)):
            if (upper_limits and upper_limits[c][i]):
                continue
            if (li_range and (bv_m[c][1][i] < li_range[0] or li_range[1] < bv_m[c][1][i])):
                continue
            if (scale_by_scatter):
                arr.append(resid[i]/scale_by_scatter(bv_m[c][1][i]))
            else:
                arr.append(resid[i])
        allClusters.append(arr)

    residual_arr = np.concatenate(allClusters)
    return np.std(residual_arr)

#takes in sigs containing two constants: constant log astrophysical scatter and constant linear measurement error
#sigs = [astrophysical scatter in logEW, measuement error in EW] ~ [.15,15]
# EW is an np array of EW values or a single number
def two_scatter(sigs,EW):
    sigma_measurement_err = np.log10(np.power(10,EW)+sigs[1]) - EW
    return np.sqrt(sigs[0]**2 + sigma_measurement_err**2)

# returns a fit to scatter as a function of equivalent width
def fit_two_scatters(bv_m,fits,upper_lim=None,omit_cluster=None):
    params = [.15,15]
    argv = (bv_m,fits,upper_lim,omit_cluster)
    res = minimize(scatter_minimizer,params,args=argv, method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing scatter."
        print "Guess: ",guess, " Result: ", res.x
    #else:
    #    print "Successful Scatter: ", res.x
    def fit(EW):
        return two_scatter(res.x,EW)
    return fit

def scatter_minimizer(params,bv_m,fits,upper_lim,omit_cluster):
    def fit(EW):
        return two_scatter(params,EW)
    return np.abs(total_scatter(bv_m,fits,omit_cluster,upper_lim,scale_by_scatter=fit) - 1)


def params_to_xy(params,X,Y,n_pin,s):
    num_x = s - 1 
    x_temp = np.concatenate((X[:max(1,n_pin)],params[:num_x],X[-1:]))
    y_temp = np.concatenate((Y[:n_pin],params[num_x:]))
    if len(x_temp) != len(y_temp):
        print s, n_pin
        print params,num_x
        print x_temp
        print y_temp
        assert len(x_temp) == len(y_temp), "Error parsing params"
    return x_temp,y_temp


# segments does not include the n_pin points
# n_pin makes the first n_pin points interpolated between
#x_method: determines how x_coordinates are selected.
#   'even' means evenly spaced in the x dimension from min(x) to max(x)
#   'bin' means spaced evenly with number per bin at least 'bins' (the next parameter) in size
#   'free' means the x_coordinates are free parameters with only the min and max fixed by the range of x 
#returns piecewise function
def general_piecewise(X,Y,segments=3,guess_fit=None, sigma=None,x_method='even', min_bin_size=2,min_length=0,monotonic=None,n_pin = 0):
    if (segments == 0): #uses a constant fit
        return np.poly1d(np.polyfit(X,Y,segments))
    if (guess_fit is None):
        guess_fit = linear_fit(X,Y)
   
    #if n_pin is 1 still starts at first point
    # x_coords defines x-locs of segment junctions starting from last pinned or start
    x_coords = np.linspace(X[max(0,n_pin-1)],X[-1],segments+1) 
    x_param = x_coords[1:-1]  # these are the free x-coordinates of the junctions
    #y_param = guess_fit(x_coords[n_pin:])
    y_param = guess_fit(x_coords) if n_pin == 0 else guess_fit(x_coords[1:])
    params0 = np.concatenate((x_param,y_param))

    def minimizer(params):
        x_temp,y_temp = params_to_xy(params,X,Y,n_pin,segments)
        #plt.plot(np.power(10,x_temp),y_temp)
        if not all(a + min_length <= b for a,b in zip(x_temp[:-1], x_temp[1:])): 
            return np.inf
        #print zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])
        if (monotonic == -1):#decreasing
            if not all(a >= b for a,b in zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])): 
                return np.inf
        if (monotonic == 1):#increasing
            if not all(a <= b for a,b in zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])):
                #return 100000 * np.sum(np.invert([a <= b for a,b in zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])]))
                return np.inf
        y_model = piecewise(x_temp,y_temp)(X)
        if sigma is None:
            return prob.chi_sqr(Y,y_model,total=True)
        return prob.chi_sqr(Y,y_model,sigma,total=True)
    
    res = minimize(minimizer,params0, method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing of %d segment piecewise function" % segments
        print "Params:",res.x
        if segments - 1 >= 0:
            print "Trying %d segment fit instead" % (segments - 1)
            return general_piecewise(X,Y,segments - 1,guess_fit,sigma,x_method, min_bin_size,min_length,monotonic,n_pin)
    
    #x_val,y_val = params_to_xy(params0,X,Y,n_pin,segments)
    x_val,y_val = params_to_xy(res.x,X,Y,n_pin,segments)
    return piecewise(x_val,y_val)

   

def pwise_fit(x,y,segments=-1,upper_lim=None, guess_fit=None, x_method='even', bin_size=MIN_PER_BIN,min_length=MIN_LENGTH):
    if (not guess_fit):
        guess_fit = linear_fit(x,y)
    if (segments == 0): #uses a constant fit
        return constant_fit(y)
    if (upper_lim == None):
        upper_lim = [False]*len(x)

    x_coords = []
    params = []
    if (x_method == 'even'):
        if (segments == -1):
            segments = 3
            x_coords = x_even(x,segments)
    elif (x_method == 'free'):
        if (segments == -1):
            x_coords = x_coords_smart_binning(x,bin_size,min_length)   
        else:
            x_coords = x_even(x,segments)
        params += x_coords[1:-1]
    elif (x_method[0] == 'b'):
        x_coords = x_coords_smart_binning(x,bin_size,min_length)
        
    sig_guess = getScatterGuess(x,y,x_coords,guess_fit)
    params += guess_fit(x_coords).tolist() + sig_guess
    argv = (x_coords,x,y,upper_lim)
    
    return combined_fit(params,argv)

def combined_fit(params,argv):
    #comp = []
    #comp.append(compare_fits(np.array(params),argv))
    res = minimize(piecewise_minimizer,params,args=argv, method='Nelder-Mead')
        
    if (not res.success):
        print "Unsuccessful minimizing of %d segment piecewise function, check initial guess" % (len(argv[0])-1)
        print res.x
    
    #comp.append(compare_fits(res.x,argv))
    #plot_comp(comp,argv[3],res.success)

    x_val,y_val,sig = get_x_y_sig(res.x,argv[0])
    #def fit(x): #returns inner function
    #        return piecewise(x,x_val,y_val)
    def scatter_fit(x):
            return step(x,x_val[1:-1],sig)
    return [piecewise(x_val,y_val), scatter_fit] #,res.success

def plot_comp(comp,upper_lim, converged):
    #plt.title("Convergence comparison")
    #plt.axhline(y=0,linestyle='--',linewidth=1,color='k')
    x = range(len(comp[0]))
    l = 'converged' if converged else 'Not converged'
    m = DOWN_ARROW if upper_lim[0] else 'o'
    plt.scatter(x[0],comp[0][0],marker=m, color='C0',label='guess fit')
    plt.scatter(x[0],comp[1][0],marker=m, color='C1',label= l)
    for i in range(1,len(x)):
        m = DOWN_ARROW if upper_lim[i] else 'o'
        plt.scatter(x[i],comp[0][i],marker=m, color='C0')
        plt.scatter(x[i],comp[1][i],marker=m, color='C1')
    
    plt.legend()
    #plt.show()
"""
def compare_fits(params,argv): #c,r,upper_lim):
    x_coords,c,r,upper_lim = argv[0],argv[1],argv[2],argv[3]
    x,y,sig = get_x_y_sig(params,x_coords)
    r_model = piecewise(x,y)(c)
    sig_model = step(c,x[1:-1],sig)
    l = []
    for i in range(len(r)):
        if (upper_lim[i]):
                #continue
                r_range = np.linspace(r_model[i] - 5*sig_model[i],r[i],100)
                l.append(np.log(np.trapz(prob.gaussian(r_model[i],sig_model[i],r_range),r_range)))
        else:
                l.append(prob.log_gaussian(r_model[i],sig_model[i],r[i]))
    return l
"""
def piecewise_minimizer(params,x_coords,c,r,upper_lim):
    x,y,sig = get_x_y_sig(params,x_coords)
    if (utils.negative_sig(sig)):
        return np.inf
    if (len(x) != len(y) or len(sig) != len(x) - 1):
        print params
        print x,y,sig
        raise RuntimeError("Error in inverse_log_likelihood")
    r_model = piecewise(x,y)(c)
    sig_model = step(c,x[1:-1],sig)
    return inverse_log_likelihood(r_model,sig_model,r,upper_lim)

#minimizing the negative of the log_likelihood
#accounts for upper limits
def inverse_log_likelihood(r_model,sig_model,r,upper_lim):
    total_sum = 0
    for i in range(len(r)):
        total_sum += norm.logpdf(r[i],loc=r_model[i],scale=sig_model[i])   ########### CHECK
        continue
        
        if (upper_lim[i]):
            total_sum += norm.logcdf(r[i],loc=r_model[i],scale=sig_model[i])
        else:
            total_sum += norm.logpdf(r[i],loc=r_model[i],scale=sig_model[i])
    return -total_sum


#takes an array of tuples (b-v, r'hk) and an operation and returns an array of tuples
# where bins of size >= bin_size tuples have been collapsed as follows: (mean b-v, op(r'hk))
def chunk(arr,bin_size,op):
        new_arr = []
        for i in range(0,len(arr),bin_size):
                if (i + 2*bin_size > len(arr)):
                        new_arr.append((np.mean([b for (b,c) in arr[i:]]), \
                                        op([c for (b,c) in arr[i:]])))
                        break
                else :
                        new_arr.append((np.mean([b for (b,c) in arr[i:i + bin_size]]), \
                                        op([c for (b,c) in arr[i:i + bin_size]])))
        return new_arr

def getScatter(c,r):
        a = copy.deepcopy(c)
        b = copy.deepcopy(r)
        detrend(a,b)
        combined = zip(a,b)
        combined.sort()
        chunked = chunk(combined,BIN_SIZE)
        [x,y] = zip(*chunked)
        return [x,y]

# finds and returns the fitline to a graph of scatter as a function of b-v color
# given an array of b-v called c and an array of r'hk called r
def getScatterFit(c,r,fit):
        a = copy.deepcopy(c)
        b = copy.deepcopy(r)
        detrend(a,b,fit)
        combined = zip(a,b)
        combined.sort()
        chunked = chunk(combined,BIN_SIZE,np.std)
        [x,y] = zip(*chunked)
        return np.poly1d(np.polyfit(x,y,1))  #what kind of fit line

def getScatterGuess(c,r,x_coords,fit):
    a = copy.deepcopy(c)
    b = copy.deepcopy(r)
    detrend(a,b,fit)
    combined = zip(a,b)
    combined.sort()
    x,y = zip(*combined)
    std = []
    for i in range(len(x_coords)-1):
        index = bisect.bisect_left(x,x_coords[i])   
        index2 = bisect.bisect_left(x,x_coords[i+1])    
        std.append(np.std(y[index:index2]))
    return std

#defines piecewise function that takes in an array of values x and fixed parameters of fit line
#returns an array for the y_values at the locations specified by x
#x_locs is increasing
def piecewise(x_locs,y_locs):
        return interpolate.interp1d(x_locs,y_locs, fill_value='extrapolate')

#x_locs defines the discontinuities,y_locs defines heights of step function
# length of x_locs is 1 fewer than y_locs
def step(x,x_locs,y_locs):
    if (type(x) != list):
        if (type(x) == float or type(x) == int):
            x = [x]
        x = x.tolist()
    res = []
    for val in x:
        i = bisect.bisect_left(x_locs,val)
        res.append(y_locs[i])
    if (len(res) == 1):
        return res[0]
    return np.array(res)

def x_coords_from_binning(c,bin_size):
    a = copy.deepcopy(c)
    a.sort()
    combined = zip(a,a)
    chunked = chunk(combined,bin_size,min)
    x,y = zip(*chunked)
    return list(y) + [a[-1]]
    
#enforces bins to have no fewer than bin_size elements and
# have length no smaller than MIN_LENGTH
def x_coords_smart_binning(c,bin_size,min_length):
    a = copy.deepcopy(c)
    a.sort()
    x_coords = [a[0]]
    i = bin_size
    while (i < len(a)):
        new_x = a[i]#np.mean(a[i-1:i+1]) #average of two adjacent star's x_coords
        if (new_x - x_coords[-1] < min_length):
            i += 1
            continue
        x_coords.append(new_x)
        i += bin_size 
    if (i == len(a)):
        x_coords.append(a[-1])
    else:
        x_coords[-1] = a[-1]
    return x_coords


# takes in the B-V values and number of segments and returns
# an evenly spaced list from min(c) to max(c)
def x_even(c,num_segments):
        return np.linspace(np.min(c),np.max(c),num_segments + 1).tolist()

"""
# chi-squared of the piecewise function which is to be minimized
#coords is [x1,x2,x3,y1,y2,y3]
#y_vals are the y values of the points specifying the piecewise fcn whose coordinates are given by x_coords
#y_coords are the free parameters changed to minimze chi2
# c in the bv color ard r is r'hk value
def chi2_piecewise(params,min_max,c,r):
        print params
        x,y,sig = get_x_y_sig(params,min_max)
        if (len(x) != len(y) or len(sig) != len(x) - 1):
                print params
                print x,y,sig
                raise RuntimeError("Error in chi2_piecewise")
        r_model = piecewise(c,x,y)
        sig_model = step(c,x[1:-1],sig)
        return np.sum(np.power(r - r_model,2)/sig_model)
"""


#if s is number of segments, there are (s-1) x_coords, (s+1) y_coords, s sigmas
#determines in x coordinates are free parameters or not
def get_x_y_sig(params,x_coords):
    params = params.tolist()
    s = len(x_coords) - 1
    if (s == len(params)/3):
            x = [x_coords[0]] + params[0:s - 1] + [x_coords[-1]]
            y = params[s - 1:2*s]
            sig = params[2*s:]
            return x,y,sig
    y = params[0:s+1]
    sig = params[s+1:]
    return x_coords,y,sig

#modifies x,y input parameters by subtracting off a linear trend
def detrend(x,y,trend):
    for i in range(len(x)):
            y[i] -= trend(x[i])

def residuals(x,y,trend):
        b = copy.deepcopy(y)
        detrend(x,b,trend)
        return b

#creates a function that interpolates between input and output.  
#returns a function
def magic_table_convert(in_column,out_column):
    if (type(in_column) == str):
        if (in_column.lower() == "teff"):
            in_column = 1
        elif (in_column.lower() == "bv" or in_column.lower() == "b-v"):
            in_column = 6
    if (type(out_column) == str):
        if (out_column.lower() == "teff"):
            out_column = 1
        elif (out_column.lower() == "bv" or out_column.lower() == "b-v"):
            out_column = 6
    
    t = ascii.read('data/mamajek_magic_table.txt')
    x,y = [],[] #x is input and y is output like a function
    for row in t:
        try:
            a = float(row[in_column])
            b = float(row[out_column])
            x.append(a)
            y.append(b)
        except:
            continue
    
    return interpolate.interp1d(x,y, fill_value='extrapolate') 


def spt_bv():
    t = ascii.read('data/mamajek_magic_table.txt')
    letters = ['O','B','A','F','G','K','M']
    #val = np.linspace(10,70,len(letters))
    #val[-1] -= 2
    spt = [row[0][:-1] for row in t[0:90]]
    sp = []
    bv = []
    for i in range(len(spt)):
    #    if (spt[i][0:2] == 'K8' or spt[i][0:2] == 'K9'):
    #        continue
        sp.append(float(spt[i][1:]) + 10*letters.index(spt[i][0]))
        bv.append(float(t[i][6]))
    
    interp = interpolate.interp1d(sp,bv, fill_value='extrapolate')
    
    #plt.plot(bv,sp)
    
    def f(x):
        return interp(utils.float_sptype(x))
    return f








#teff is an array
#At every T in teff it uses the soderblom 1993 Pleiades table to find the log(EW) corresponding to
# PRIMORDIAL_NLI = 3.2
# returns an array
def teff_to_primli(teff):
    t = genfromtxt('data/NLi_to_LiEW.csv', delimiter=',')
    logEW = [row[0] for row in t[1:]]
    #print logEW
    arr = []
    for temp in teff:
        #make my own column via interp
        col = []
        for row in t[1:]:
            #print row[1:]
            #print t[0][1:]
            col.append(interpolate.interp1d(t[0][1:], row[1:],fill_value='extrapolate')(temp))
        arr.append(interpolate.interp1d(col,logEW, fill_value='extrapolate')(PRIMORDIAL_NLI).tolist())
            
    return arr

#teff is an array
#At every T,nli it uses the soderblom 1993 Pleiades table to find the log(EW)
# returns an array
def teff_nli_to_li(teff,NLI):
    t = genfromtxt('data/NLi_to_LiEW.csv', delimiter=',')
    t2 = genfromtxt('data/zapatero_osorio_teff_nli_ewli.txt')
    t2[1:,1:] = np.log10(1000*t2[1:,1:]) #convert to logEW
    logEW = [row[0] for row in t[1:]]
    temp_axis = [row[0] for row in t2[1:]]
    #print logEW
    arr = []
    for temp,nli in zip(teff,NLI):
        #make my own column via interp of exact temp
        col = []
        if temp <= 4000:
            for row in t2[1:]:
                #print row[1:]
                #print t2[0][1:]
                col.append(interpolate.interp1d(t2[0][1:], row[1:],fill_value='extrapolate')(nli))
            arr.append(interpolate.interp1d(temp_axis,col,fill_value='extrapolate')(temp).tolist())
        else:
            for row in t[1:]:
                #print row[1:]
                #print t[0][1:]
                col.append(interpolate.interp1d(t[0][1:], row[1:],fill_value='extrapolate')(temp))
            arr.append(interpolate.interp1d(col,logEW, fill_value='extrapolate')(nli).tolist())
    return np.array(arr)

#return array of primordial Li for each B-V value in li_constants.BV
def primordial_li(ngc2264_fit=None,fromFile=True, saveToFile=False):
    assert (ngc2264_fit or fromFile),"primordial_li must take in ngc2264 fit if not reading from a file"
    import li_constants as const
    if (fromFile):
       prim_li = pickle.load(open('data/primordial_li.p','rb'))
       return interpolate.interp1d(const.BV,prim_li, fill_value='extrapolate')
    
    teff = magic_table_convert(6,1)(const.BV) #convert B-V to Teff
    prim_li = teff_to_primli(teff)
    
    SODERBLOM_4000K_BOUNDARY = 1.356 #B-V
    BTNEXTGEN_BOUNDARY = 1.522 # No lithium has decreased by 5 Myr
    loc1 = bisect.bisect_left(const.BV,SODERBLOM_4000K_BOUNDARY)
    loc2 = bisect.bisect_left(const.BV,BTNEXTGEN_BOUNDARY)
    #print [prim_li[loc1],ngc2264_fit(BTNEXTGEN_BOUNDARY)]
    middle = interpolate.interp1d([SODERBLOM_4000K_BOUNDARY,BTNEXTGEN_BOUNDARY],[prim_li[loc1],ngc2264_fit(BTNEXTGEN_BOUNDARY)], fill_value='extrapolate')(const.BV[loc1:loc2])

    final_li = prim_li[0:loc1] + middle.tolist() + ngc2264_fit(const.BV[loc2:]).tolist()

    #print final_li
    if (saveToFile):
        pickle.dump(final_li,open('data/primordial_li.p','wb'))
    return interpolate.interp1d(const.BV,final_li, fill_value='extrapolate')

# takes in the logR'HK value and returns the age in units Myr
#def getMamaAge(r):
#    return np.power(10,-38.053 - 17.912*r - 1.6675*r*r)/1e6

# Takes in age in units Myr and converts to logR'HK
#def getMamaRHK(age):
#    log_age = np.log10(np.array(age)*1e6)
#    return 8.94 - 4.849*log_age + .624*log_age**2 - .028 * log_age**3

def VI_to_teff(in_column=2,out_column=6,switch=False):
    #if switch:
    #    in_column,out_column = out_column,in_column
    t = ascii.read('data/Teff_to_V-I.csv')
    x,y = [],[] #x is input and y is output like a function
    for row in t:
        try:
            a = float(row[in_column])
            b = float(row[out_column])
            if np.isnan(a) or np.isnan(b):
                continue
            x.append(a)
            y.append(b)
        except:
            continue
    VI = interpolate.interp1d(x,y, fill_value='extrapolate') 
    VI2 = lambda x: 9581.1 + -14264*x + 40759*x**2 - 74141*x**3 + 60932*x**4 - 18021*x**5
    return lambda vi: (vi < 1.2) * VI2(vi) + (vi >= 1.2) * VI(vi)

# updates fits
def cluster_scatter_from_stars(bv_m,fits):
    const = utils.init_constants('lithium')
    bv_threshold = 0.03
    arr = []
    for c in range(len(fits)):
        num_stars = []
        bv_arr = np.array(bv_m[c][0])
        for bv in const.BV:
            mask = (bv_arr >= bv - bv_threshold) & (bv_arr <= bv + bv_threshold)
            num_stars.append(np.sum(mask))
        arr.append(num_stars)
    arr = np.array(arr)
    
    arr = arr / (np.sum(arr,axis=0) + .001)
    #arr holds % of stars each cluster has at each B-V slice
    
    arr = (1 - arr)*.35 + 0.05

    for c in range(len(arr)):
        max_bv,min_bv = max(bv_m[c][0]),min(bv_m[c][0])
        dist = np.minimum(np.abs(const.BV-max_bv),np.abs(const.BV-min_bv))
        dist = dist * ((const.BV > max_bv) | (const.BV < min_bv))
        range_offset = (0.4/0.2)*dist
        arr[c] += range_offset

    arr = savgol_filter(arr, 51, 3)
    
    for c in range(len(fits)):
        fits[c][1] = piecewise(const.BV,arr[c])
    return fits



def get_fit_residuals(bv_m,fits,metal,upper_limits=None,li_range=None,age_range=None,linSpace=False,scale_by_std=False):
    const = utils.init_constants(metal)
    allClusters = []
    #residual_arr = []

    for c in range(len(fits)):
        if age_range is not None and not (age_range[0] <= const.CLUSTER_AGES[c]\
                <= age_range[1]):
            allClusters.append([])
            continue
        arr = []
        resid = None
        if linSpace:
            resid = np.power(10,bv_m[c][1]) - np.power(10,fits[c][0](bv_m[c][0]))
        else: #log space
            resid = residuals(bv_m[c][0],bv_m[c][1],fits[c][0])  #Log space
        for i in range(len(resid)):
            if (upper_limits is not None and upper_limits[c][i]):
                continue
            if (li_range is not None and (bv_m[c][1][i] < li_range[0] or \
                    li_range[1] < bv_m[c][1][i])):
                continue
            arr.append(resid[i])
        if scale_by_std:
            arr = np.array(arr)/np.std(arr)
        allClusters.append(arr)

    residual_arr = np.concatenate(allClusters)
    return allClusters,residual_arr

def fit_histogram(metal,residual_arr=None,fromFile=True,saveToFile=False):
    if fromFile:
        [x,pdf,cdf] = np.load('grids/' + metal + '_likelihood_fit.npy')
        return piecewise(x,pdf),piecewise(x,cdf)
    const = utils.init_constants(metal)
    
    assert residual_arr is not None or fromFile, "Must provide residuals if not \
            reading from file"
    mu = np.mean(residual_arr)
    sigma = np.std(residual_arr)

    x = np.linspace(np.min(residual_arr)-.5,np.max(residual_arr)+.5,1000) #1000 for linear?
    lim = 2
    if metal=='calcium':
        lim = 5
        x = np.linspace(np.min(residual_arr)-.5,np.max(residual_arr)+.1,800)
    before,after = np.linspace(-lim,min(x),50),np.linspace(max(x),lim,50)
    x = np.concatenate((before,x,after))
    cdf = np.array([(residual_arr < n).sum() for n in x],dtype='float')/len(x)
    cdf /= cdf[-1]
    #plt.plot(x,cdf,label='cdf_og')

    if metal=='calcium':
        smoothed = savgol_filter(cdf, 55, 3)
        smoothed = savgol_filter(smoothed, 25, 3)
        smoothed = savgol_filter(smoothed, 9, 3)
    else:
        smoothed = savgol_filter(cdf, 55, 3)
        smoothed = savgol_filter(smoothed, 25, 3)
        smoothed = savgol_filter(smoothed, 9, 3)
   
    #plt.plot(x,smoothed,label='smoothed')
    #plt.plot(x,norm.cdf(x,loc=mu,scale=sigma),label='gaussian cdf')
    #plt.legend()
    #plt.show()

    pdf = np.gradient(smoothed)
    prob.normalize(x,pdf)
    
    #plt.plot(x,pdf,label='smoothed')
    #plt.show()
    
    inds = np.nonzero(pdf > 1.5)[0] if metal=='lithium' else \
            np.nonzero(pdf > max(pdf)/2)[0]
    i,j = inds[0],inds[-1]
    def exp_fit(x,a,b,c):
        return a*np.exp(b*x + c)
    popt,pcov = curve_fit(exp_fit,x[:i],pdf[:i],p0=[5,5,-1])
    pdf[:i] = exp_fit(x[:i],*popt)
    popt,pcov = curve_fit(exp_fit,x[j:],pdf[j:],p0=[.5,-9,2.5])
    pdf[j:] = exp_fit(x[j:],*popt)
    
    #pdf[:i] = savgol_filter(pdf[0:i], 55, 3)
    #pdf[j:] = savgol_filter(pdf[j:], 55, 3)
    #pdf[i-10:i+10] = savgol_filter(pdf[i-10:i+10], 9, 3)
    #pdf[j-20:j+20] = savgol_filter(pdf[j-20:j+20], 21, 3)
    pdf = savgol_filter(pdf, 9, 3)
    
    if metal=='calcium':
        m,n = np.nonzero(pdf >= 0.32)[0][0],np.nonzero(pdf >= 0.45)[0][0]
        #plt.plot(x,pdf)
        #plt.show()
        pdf[m:n] = piecewise([x[m],x[n]],[pdf[m],pdf[n]])(x[m:n])
        
        pdf = savgol_filter(pdf,21,3)
        #plt.plot(x,pdf)
        #plt.show()
    
    pdf [:2] = [0,0]
    pdf[-2:] = [0,0]
    prob.normalize(x,pdf)
    
    cdf = integrate.cumtrapz(pdf, x=x, initial=0)
    cdf /= cdf[-1]
    
    #plt.plot(x,cdf,label='cdf')
    #plt.plot(x,norm.cdf(x,loc=mu,scale=sigma),label='gaussian cdf')
    #plt.legend()
    #plt.show()

    if saveToFile:
        np.save('grids/' + metal + '_likelihood_fit',[x,pdf,cdf])
    return piecewise(x,pdf),piecewise(x,cdf)




#######################################################
# ldb fit is the bldb fit but doesn't conflict with name
# only removes omit_cluster for lithium, vs_age_fits does it for calcium only
def get_valid_metal(bv,fits,const,primordial_li_fit=None,ldb_fit=None,omit_cluster=None):
    #if const.METAL_NAME == 'calcium': #### SHORTCUT #####
        #rhk = [fits[i][0](bv) for i in range(len(fits))]
        #sig = [fits[i][1](bv) for i in range(len(fits))] #Change if non-constant scatter
        #if omit_cluster is not None:
        #    del rhk[omit_cluster]
        #    del sig[omit_cluster]
        #    del const.CLUSTER_AGES[omit_cluster]
        #    del const.CLUSTER_NAMES[omit_cluster]
        #return rhk,sig,const.CLUSTER_AGES,const.CLUSTER_NAMES

    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = [],[],[],[]
    #info added from primordial lithium
    if (const.METAL_NAME == 'lithium'):
        if primordial_li_fit is None:
            primordial_li_fit = MIST_primordial_li()#fits[0][0],fromFile=False,saveToFile=True)
        CLUSTER_AGES.append(const.PRIMORDIAL_LI_AGE)
        rhk.append(float(primordial_li_fit(bv)))
        CLUSTER_NAMES.append("Primordial EWLi")
        scatter.append(0.05)
    for i in range(len(fits)):
        if (omit_cluster is not None and i == omit_cluster and const.METAL_NAME=='lithium'):
            continue
        r = fits[i][0](bv)
        if (const.METAL_RANGE[0] <=  r <= const.METAL_RANGE[1]):
            rhk.append(r)
            scatter.append(fits[i][1](bv))
            CLUSTER_AGES.append(const.CLUSTER_AGES[i])
            CLUSTER_NAMES.append(const.CLUSTER_NAMES[i])
    #info added from depletion boundary
    if (const.METAL_NAME == 'lithium' and bv >= const.BLDB_LOWER_LIM):
        if ldb_fit is None:
            ldb_fit = bldb_fit(fits)
        scatter.append(0.15)
        CLUSTER_NAMES.append('BLDB point')
        CLUSTER_AGES.append(ldb_fit(bv))
        rhk.append(const.ZERO_LI)
    return rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES

#def constr_inverse_x(x1,y1):
#    def inverse_x(x,a,b,c):
#        e = y1 - a/(b*x1 + c)
#        return a/(b*x + c) + e
#    return inverse_x

#def make_half_n_half(X,Y):
#    def half_n_half(x,x2):
#        y2 = poly_fit(X[:4],Y[:4],1)(x2)
#        linear_fit = poly_fit([X[0],x2],[Y[0],y2],1)
#        constr_poly_fit = constrained_poly_fit(X,Y,x2,y2)
#        return np.where(x < x2,linear_fit(x),constr_poly_fit(x))
#    return half_n_half

#def simple_pwise_fit(x,y,segments):
#    lin_fit = linear_fit(x[:2],y[:2])
#    x,y,segments = x[1:],y[1:],segments - 1
#    guess_fit = constrained_poly_fit(x,y)
#    x_coords = np.linspace(np.min(x),np.max(x),segments + 1)
#    y_coords = guess_fit(x_coords)
#    p0 = np.concatenate((x_coords[1:-1],y_coords[1:]))
#    def f(X,*params):
#        x_c = np.concatenate(([x[0]],params[:segments - 1],[x_coords[-1]]))
#        y_c = np.concatenate(([y[0]],params[segments - 1:]))
#        return piecewise(x_c,y_c)(X)
#        
#        #return piecewise(x_coords,np.insert(params,0,y[0]))(X)
#    popt,pcov = curve_fit(f,x,y,p0=p0)
#    return lambda age : np.where(np.array(age) < x[0],lin_fit(age),f(age,*popt))

#returns functions that take in age to get metal,sigma
def vs_age_fits(bv,cluster_ages,rhk,scatter,metal,li_scatter_fit=None,omit_cluster=None):
    if not li_scatter_fit and metal == 'lithium':
        assert False
        return
    cluster_ages_og,scatter_og = copy.deepcopy(cluster_ages),copy.deepcopy(scatter)
    if omit_cluster is not None and metal == 'calcium': 
        del cluster_ages[omit_cluster]
        del rhk[omit_cluster]
        del scatter[omit_cluster]
    metal_fit = None
    mu_lbl = ''
    if (metal == 'calcium'):
        const = utils.init_constants('calcium')
        num_stars = np.array(const.CLUSTER_INDEX)[:,1] - np.array(const.CLUSTER_INDEX)[:,0]
        err = 1.0/num_stars
        metal_fit = constrained_poly_fit(np.log10(cluster_ages),rhk,lim=0,sigma=err)
        #metal_fit = poly_fit(np.log10(cluster_ages),rhk,2,weight=num_stars)
        mu_lbl = '2nd-order Polynomial'
    elif (metal == 'lithium'):
        cluster_ages = np.log10(cluster_ages)
        bv_cut = [0.65,1,1.6]
        segs = [2,3,2,1]
        s = segs[bisect.bisect_left(bv_cut,bv)]
        
        metal_fit = general_piecewise(cluster_ages,rhk,s,\
                n_pin=2,monotonic=-1,min_length=.2,sigma=scatter)
            
    def mu(age):
        return metal_fit(np.log10(age))

    sig = None
    #find fit to mu of ca/li for grid
    #5 different methods for handling scatter:
    # gaussian fit,total detrended mean, mean clusters, best-fit, linear interp
    if (metal == 'calcium'):
        sig = fit_gaussian(np.log10(cluster_ages_og),scatter_og)
    elif (metal == 'lithium'):
        def scatter(age):
            return li_scatter_fit(mu(age))
        sig = scatter

    return mu,sig,mu_lbl


###########################################

#return array of primordial Li for each B-V value in li_constants.BV
def MIST_primordial_li(ngc2264_fit=None,fromFile=True, saveToFile=False):
    assert (ngc2264_fit or fromFile),"primordial_li must take in ngc2264 fit if not reading from a file"
    import li_constants as const
    if (fromFile):
       prim_li = pickle.load(open('data/mist_primordial_li.p','rb'))
       return interpolate.interp1d(const.BV,prim_li, fill_value='extrapolate')
    
    teff = magic_table_convert('bv','teff')(const.BV) #convert B-V to Teff
   
    t1 = ascii.read('data/MIST_iso_1Myr.txt')
    t5 = ascii.read('data/MIST_iso_5Myr.txt')
    
    star_mass5 = interpolate.interp1d(t5['log_Teff'][0:275],t5['initial_mass'][0:275],\
            fill_value='extrapolate')(np.log10(teff))
    Nli5 = 12 + np.log10(interpolate.interp1d(t5['log_Teff'][0:275],t5['surface_li7'][0:275],\
            fill_value='extrapolate')(np.log10(teff)))
    
    Nli1 = 12 + np.log10(interpolate.interp1d(t1['initial_mass'],t1['surface_li7'],\
            fill_value='extrapolate')(star_mass5))
    #Nli1 = 12 + np.log10(interpolate.interp1d(t1['log_Teff'],t1['surface_li7'],\
    #        fill_value='extrapolate')(np.log10(teff)))
   
    #convert teff and nli to EW   
    deltaEW = teff_nli_to_li(teff,Nli1) - teff_nli_to_li(teff,Nli5)
    deltaEW = np.clip(deltaEW,0,None)
    assert all(deltaEW >= 0), "EW must decrease monotonically"
    final_li = ngc2264_fit(const.BV) + deltaEW
    
    #final_li = teff_nli_to_li(teff,[3.2]*len(teff))
    #final_li = teff_nli_to_li(teff,Nli5)

    if (saveToFile):
        pickle.dump(final_li,open('data/mist_primordial_li.p','wb'))
    return interpolate.interp1d(const.BV,final_li, fill_value='extrapolate')







# refreshes everything saved to files
def main():
    upper_lim = None
    bv_m,fits = readData.read_calcium()
    _,res_arr = get_fit_residuals(bv_m,fits,'calcium',upper_lim,li_range=None,linSpace=False,scale_by_std=True)
    fit_histogram('calcium',residual_arr=res_arr,fromFile=False,saveToFile=True)
    
    
    bv_m, upper_lim, fits = readData.read_lithium()
    MIST_primordial_li(ngc2264_fit=fits[0][0],fromFile=False, saveToFile=True)
    _,res_arr=get_fit_residuals(bv_m,fits,'lithium',upper_lim,li_range=None,linSpace=False)
    fit_histogram('lithium',residual_arr=res_arr,fromFile=False,saveToFile=True)
    





if  __name__ == "__main__":
    main()










