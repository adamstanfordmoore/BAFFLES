"""
Adam Stanford-Moore
5/22/20
Module giving necessary fitting functions
"""

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate,integrate
from scipy.stats import t as scipy_t
from scipy.stats import norm
from scipy.optimize import minimize,curve_fit
from scipy.signal import savgol_filter
import time
import copy
import bisect
import pickle
from astropy.io import ascii
import probability as prob
import utils
import readData
import plotting as my_plot
from os.path import join

BIN_SIZE = 10
MIN_PER_BIN = 4
MIN_LENGTH = .15 #minimum B-V length of piecewise segment
#DOWN_ARROW = u'$\u2193$'
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
        print("Unsuccessful minimizing of constrained polynomial function, check initial guess")
        print(res.x)

    return poly_constr_vert(res.x,x1,y1) if lim is None else np.poly1d(res.x)

#def right_cubic_root(params):
#    root = 4*params[1]*params[1] - 12*params[0]*params[2]
#    if (root < 0):
#        return -100000
#        #raise RuntimeError("CUBIC ERROR")
#    a = (-2*params[1] - np.sqrt(root)) / (6 * params[0])
#    b = (-2*params[1] + np.sqrt(root)) / (6 * params[0])
#    if (a >= b):
#        return a
#    return b
    
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
   
    upper_clusters = ['M67','Hyades','M34','M35']
    upper_ages = [ages[cluster_names.index(name)] for name in upper_clusters]
    upper_bv = [bv_at_zero_li[cluster_names.index(name)] for name in upper_clusters]
    fit = piecewise(upper_bv,np.log10(upper_ages))
    
    def ldb_age(x):
        return np.power(10,fit(x))
    
    if (plot):
        pp = PdfPages('plots/bldb_vs_age.pdf')
        plt.xlabel('B-V',size=my_plot.AXIS_LABEL_SIZE)
        plt.ylabel('Age (Myr)',size=my_plot.AXIS_LABEL_SIZE)
        plt.yscale('log')
        for c in range(len(ages)):
            plt.scatter(bv_at_zero_li[c],ages[c],s=60,label=cluster_names[c],color=const.COLORS[c+1],marker=const.MARKERS[c+1])
        plt.plot(BV,ldb_age(BV))
        plt.fill_between(BV,ldb_age(BV),color='C0',alpha=.2,label="valid ages")
        plt.legend()        
        plt.minorticks_on()
        plt.tick_params(axis='both',which='both',right=True,top=True)
        plt.tight_layout()
        pp.savefig()
        #plt.show()
        pp.close()
        plt.close()
    
    return ldb_age

#returns the 1d polynomial
def linear_fit(x,y,scatter=False):
    return poly_fit(x,y,1,scatter=scatter)

def constant_fit(y):
    m,s = np.mean(y),np.std(y)
    return [np.poly1d([m]),np.poly1d([s])]

#def dip_gaussian(x,mu,sig,A,C):
#    chi2 = prob.chi_sqr(x,mu,sig)
#    return C + A/(sig*np.sqrt(2*np.pi)) * np.exp(-chi2/2)



# dip_bv_range and dip_li_range specify box containing dipped stars,
# so they can be excluded for background polynomial
# total_dip_fit_range specifies bv range to fit polynomial
# edge_box specifies box of stars (bv_min,bv_max,li_min,li_max) to make dip go through
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

    no_dip_fit,no_dip_sig_fit = poly_fit(bv_,li_,2,upper_lim_)
    
    edge_mask = (bv_dip > edge_box[0]) & (bv_dip < edge_box[1]) & \
            (li_dip > edge_box[2]) & (li_dip < edge_box[3])

    sigma = np.where(edge_mask,.02,1) # makes fit go through edge_box points

    def dip_poly(x,a,b,c):
        return a*x**2 + b*x + c

    popt,pcov = curve_fit(dip_poly,bv_dip,li_dip,sigma=sigma)#,p0=[100,-90,21.25])
    #print(popt)
    #popt,pcov = curve_fit(dip_gaussian,bv_dip,li_dip,p0=[0.45,.04,-.1],sigma=sigma)
    dip_fit = lambda x: dip_poly(np.array(x),*popt)
    #dip_fit = lambda x: dip_gaussian(np.array(x),*popt)
    dip_sig_fit = np.poly1d([.25])

    final_fit = lambda x : np.where((np.array(x) <= dip_bv_range[1]) & \
            (np.array(x) >= dip_bv_range[0]),dip_fit(x),no_dip_fit(x))
    
    final_sig_fit = lambda x: np.where((np.array(x) <= dip_bv_range[1]) & \
            (np.array(x) >= dip_bv_range[0]),dip_sig_fit(x),no_dip_sig_fit(x))
    
    return [final_fit,final_sig_fit]


# finds and returns n dimensional polynomial fit.
#If scatter=True than it returns [median fit,scatter fit]
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
        print("Unsuccessful minimizing of polynomial function, check initial guess")
        print(res.x)

    sig = res.x[-1]
    poly = res.x[0:-1]
    fit = np.poly1d(poly)
    return [fit, constant_fit(sig)[0]]

def poly_minimizer(params,c,r,upper_lim):
    sig = params[-1]
    if (utils.negative_sig(sig)):
        return np.inf
    r_model = np.poly1d(params[0:-1])(c)
    sig_model = constant_fit(sig)[0](c)       
    return inverse_log_likelihood(r_model,sig_model,r,upper_lim)

# 4 params = mu,sig,amplitude,vertical_offset
# takes in log ages and returns function of normal ages
def fit_gaussian(x,y):
    guess = [np.mean(x),np.std(x),np.max(y)-np.min(y),np.min(y)]
    res = minimize(gaussian_scatter_minimizer,guess,args=(x,y), method='Nelder-Mead')
    if (not res.success):
        print("Unsuccessful minimizing of gaussian scatter_vs_age fit")
        print(res.x)

    [mu,sig,A,c] = res.x
    #print("Gaussian fit params= ", res.x)
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
    sigma_weight = 1.0/num_stars
    return prob.chi_sqr(y,fit,sigma_weight,total=True)

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

    def params_to_xy(params,X,Y,n_pin,s):
        num_x = s - 1
        x_temp = np.concatenate((X[:max(1,n_pin)],params[:num_x],X[-1:]))
        y_temp = np.concatenate((Y[:n_pin],params[num_x:]))
        if len(x_temp) != len(y_temp):
            print(s, n_pin)
            print(params,num_x)
            print(x_temp)
            print(y_temp)
            assert len(x_temp) == len(y_temp), "Error parsing params"
        return x_temp,y_temp

    def minimizer(params):
        x_temp,y_temp = params_to_xy(params,X,Y,n_pin,segments)
        if not all(a + min_length <= b for a,b in zip(x_temp[:-1], x_temp[1:])): 
            return np.inf
        if (monotonic == -1):#decreasing
            if not all(a >= b for a,b in zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])): 
                return np.inf
        if (monotonic == 1):#increasing
            if not all(a <= b for a,b in zip(y_temp[max(0,n_pin-1):-1], y_temp[max(1,n_pin):])):
                return np.inf
        y_model = piecewise(x_temp,y_temp)(X)
        if sigma is None:
            return prob.chi_sqr(Y,y_model,total=True)
        return prob.chi_sqr(Y,y_model,sigma,total=True)
    
    res = minimize(minimizer,params0, method='Nelder-Mead')
    if (not res.success):
        print("Unsuccessful minimizing of %d segment piecewise function" % segments)
        print("Params:",res.x)
        if segments - 1 >= 0:
            print("Trying %d segment fit instead" % (segments - 1))
            return general_piecewise(X,Y,segments - 1,guess_fit,sigma,x_method, min_bin_size, \
                                     min_length,monotonic,n_pin)
    
    x_val,y_val = params_to_xy(res.x,X,Y,n_pin,segments)
    return piecewise(x_val,y_val)

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
    spt = [row[0][:-1] for row in t[0:90]]
    sp = []
    bv = []
    for i in range(len(spt)):
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
    arr = []
    for temp in teff:
        #make my own column via interp
        col = []
        for row in t[1:]:
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
    #print(logEW
    arr = []
    for temp,nli in zip(teff,NLI):
        #make my own column via interp of exact temp
        col = []
        if temp <= 4000:
            for row in t2[1:]:
                col.append(interpolate.interp1d(t2[0][1:], row[1:],fill_value='extrapolate')(nli))
            arr.append(interpolate.interp1d(temp_axis,col,fill_value='extrapolate')(temp).tolist())
        else:
            for row in t[1:]:
                col.append(interpolate.interp1d(t[0][1:], row[1:],fill_value='extrapolate')(temp))
            arr.append(interpolate.interp1d(col,logEW, fill_value='extrapolate')(nli).tolist())
    return np.array(arr)

def VI_to_teff(in_column=2,out_column=6):
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

def get_fit_residuals(bv_m,fits,metal,upper_limits=None,li_range=None,age_range=None,
                      linSpace=False,scale_by_std=False,vs_age_fit=True,zero_center=True):
    const = utils.init_constants(metal)
    allClusters = []

    grid_median = np.load(const.DEFAULT_MEDIAN_GRID)
    mu_interp = interpolate.interp2d(const.AGE,const.BV_S,grid_median) if metal == 'lithium' else \
                      interpolate.interp1d(const.AGE,grid_median)

    for c in range(len(fits)):
        if age_range is not None and not (age_range[0] <= const.CLUSTER_AGES[c]\
                <= age_range[1]):
            allClusters.append([])
            continue
        arr = [] #holds non UL from cluster i
        resid = None
        if vs_age_fit:
            resid = None
            if metal == 'lithium':
                resid = np.array(bv_m[c][1]) - mu_interp(const.CLUSTER_AGES[c],bv_m[c][0]).flatten()
            else:        
                resid = np.array(bv_m[c][1]) - mu_interp(const.CLUSTER_AGES[c])
        elif linSpace:
            resid = np.power(10,bv_m[c][1]) - np.power(10,fits[c][0](bv_m[c][0]))
        else: #log space
            resid = residuals(bv_m[c][0],bv_m[c][1],fits[c][0])  #Log space
        
        #now filter out upper-limits from resid
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
    if zero_center:
        #print("Subtracting off median of %.5f from residuals, (mean = %.5f)" % (np.median(residual_arr),np.mean(residual_arr)))
        resid_mean = np.median(residual_arr) 
        residual_arr -= resid_mean
        for i in range(len(allClusters)):
            allClusters[i] = np.array(allClusters[i]) - resid_mean
    return allClusters,residual_arr


def lorentz_pdf(x,xo,gamma):
    return gamma/((x-xo)**2 + gamma**2)/np.pi

def lorentz_cdf(x,xo,gamma):
    return np.arctan((x-xo)/gamma)/np.pi + 0.5
 
def student_pdf(x,df, scale):
    return scipy_t.pdf(x, df, 0, scale)
   
def student_cdf(x,df, scale):
    return scipy_t.cdf(x, df, 0, scale)


def fit_student_t(metal,residual_arr=None,fromFile=True,saveToFile=False):
    assert residual_arr is not None or fromFile, "Must provide residuals if not \
            reading from file"
    popt = None
    if fromFile:
        popt = np.load(join('grids/',metal + '_student_t_likelihood_params.npy'))
    else:
        buf = 0.1
        x = np.linspace(np.min(residual_arr)-buf,np.max(residual_arr)+buf,1000) #1000 for linear?
        cdf = np.array([(residual_arr < n).sum() for n in x],dtype='float')
        cdf /= cdf[-1]
    
        #popt,pcov = curve_fit(lorentz_cdf,x,cdf)
        popt,pcov = curve_fit(student_cdf,x,cdf,p0=[1,.1])
    
        print("Optimal params=", popt)  
    
    def pdf_fit(input):
        return student_pdf(input, *popt) 
        #return lorentz_pdf(input,*popt)
    def cdf_fit(input):
        return student_cdf(input, *popt) 
        #return lorentz_cdf(input,*popt)

    if not fromFile and saveToFile:
        np.save(join('grids', metal + '_student_t_likelihood_params'),popt)
    return pdf_fit,cdf_fit


def fit_histogram(metal,residual_arr=None,fromFile=True,saveToFile=False):
    if fromFile:
        [x,pdf,cdf] = np.load(join('grids',metal + '_likelihood_fit.npy'))
        return piecewise(x,pdf),piecewise(x,cdf)
    const = utils.init_constants(metal)
    
    assert residual_arr is not None or fromFile, "Must provide residuals if not \
            reading from file"
    mu = np.mean(residual_arr)
    sigma = np.std(residual_arr)

    x = np.linspace(np.min(residual_arr)-.5,np.max(residual_arr)+.5,1000) #1000 for linear?
    lim = 2
    if metal=='calcium':
        #lim = 5
        lim = 1
        x = np.linspace(np.min(residual_arr)-.5,np.max(residual_arr)+.1,800)
    before,after = np.linspace(-lim,min(x),50),np.linspace(max(x),lim,50)
    x = np.concatenate((before,x,after))
    cdf = np.array([(residual_arr < n).sum() for n in x],dtype='float')
    cdf /= cdf[-1]

    if metal=='calcium':
        smoothed = savgol_filter(cdf, 145, 3)
        smoothed = savgol_filter(smoothed, 55, 3)
    else:
        smoothed = savgol_filter(cdf, 85, 3)
        smoothed = savgol_filter(smoothed, 55, 3)
        #smoothed = savgol_filter(cdf, 55, 3)
        #smoothed = savgol_filter(smoothed, 25, 3)
        #smoothed = savgol_filter(smoothed, 9, 3)
   
    pdf = np.gradient(smoothed)
    prob.normalize(x,pdf)

    inds = np.nonzero(pdf > max(pdf)/2)[0]
    #inds = np.nonzero(pdf > 1.5)[0] if metal=='lithium' else \
    #        np.nonzero(pdf > max(pdf)/2)[0]
    i,j = inds[0],inds[-1]
    def exp_fit(x,a,b,c):
        return a*np.exp(b*x + c)
    popt,pcov = curve_fit(exp_fit,x[:i],pdf[:i],p0=[5,5,-1])
    pdf[:i] = exp_fit(x[:i],*popt)
    popt,pcov = curve_fit(exp_fit,x[j:],pdf[j:],p0=[.5,-9,2.5])
    pdf[j:] = exp_fit(x[j:],*popt)
    pdf = savgol_filter(pdf, 9, 3)

    if metal=='calcium':
        m,n = np.nonzero(pdf >= 0.32)[0][0],np.nonzero(pdf >= 0.45)[0][0]
        pdf[m:n] = piecewise([x[m],x[n]],[pdf[m],pdf[n]])(x[m:n])        
        pdf = savgol_filter(pdf,21,3)
    
    pdf [:2] = [0,0]
    pdf[-2:] = [0,0]
    prob.normalize(x,pdf)    
    cdf = integrate.cumtrapz(pdf, x=x, initial=0)
    cdf /= cdf[-1]

    if saveToFile:
        np.save(join('grids',metal + '_likelihood_fit'),[x,pdf,cdf])
    return piecewise(x,pdf),piecewise(x,cdf)




#######################################################
# ldb fit is the bldb fit but doesn't conflict with name
# only removes omit_cluster for lithium, vs_age_fits does it for calcium only
def get_valid_metal(bv,fits,const,primordial_li_fit=None,ldb_fit=None,omit_cluster=None):
    rhk,scatter,CLUSTER_AGES,CLUSTER_NAMES = [],[],[],[]
    #info added from primordial lithium
    if (const.METAL_NAME == 'lithium'):
        if primordial_li_fit is None:
            primordial_li_fit = MIST_primordial_li()#fits[0][0],fromFile=False,saveToFile=True)
        CLUSTER_AGES.append(const.PRIMORDIAL_LI_AGE)
        rhk.append(float(primordial_li_fit(bv)))
        CLUSTER_NAMES.append("Primordial LiEW")
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

#returns functions that take in age to get metal,sigma
def vs_age_fits(bv,cluster_ages,rhk,scatter,metal,omit_cluster=None):
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
        if omit_cluster is not None:
            err = np.delete(err,omit_cluster)
        metal_fit = constrained_poly_fit(np.log10(cluster_ages),rhk,lim=0,sigma=err)
        #metal_fit = poly_fit(np.log10(cluster_ages),rhk,2,weight=num_stars)
        mu_lbl = 'polynomial fit'
    elif (metal == 'lithium'):
        cluster_ages = np.log10(cluster_ages)
        bv_cut = [0.65,1,1.6]
        segs = [2,3,2,1]
        s = segs[bisect.bisect_left(bv_cut,bv)]
        
        # hardcode fit to go through hyades during the lithium fit
        if 0.41 < bv < 0.51:
            const = utils.init_constants('lithium')
            scatter[const.CLUSTER_NAMES.index("Hyades") + 1] = 0.01
            scatter = scatter[:-1] #del scatter[const.CLUSTER_NAMES.index("M67") + 1]
            cluster_ages = cluster_ages[:-1]#[const.CLUSTER_NAMES.index("M67") + 1]
            rhk = rhk[:-1] # [const.CLUSTER_NAMES.index("M67") + 1]
            #s += 1

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
        sig = np.poly1d([np.mean(scatter)]) #placeholder--isn't used

    return mu,sig,mu_lbl


###########################################

#return array of primordial Li for each B-V value in li_constants.BV
def MIST_primordial_li(ngc2264_fit=None,fromFile=True, saveToFile=False):
    assert (ngc2264_fit or fromFile),"primordial_li must take in ngc2264 fit if not reading from a file"
    import li_constants as const
    if (fromFile):
       prim_li = pickle.load(open(join('data','mist_primordial_li.p'),'rb'))
       return interpolate.interp1d(const.BV,prim_li, fill_value='extrapolate')
    
    teff = magic_table_convert('bv','teff')(const.BV) #convert B-V to Teff
   
    t1 = ascii.read(join('data','MIST_iso_1Myr.txt'))
    t5 = ascii.read(join('data','MIST_iso_5Myr.txt'))
    
    star_mass5 = interpolate.interp1d(t5['log_Teff'][0:275],t5['initial_mass'][0:275],\
            fill_value='extrapolate')(np.log10(teff))
    Nli5 = 12 + np.log10(interpolate.interp1d(t5['log_Teff'][0:275],t5['surface_li7'][0:275],\
            fill_value='extrapolate')(np.log10(teff)))
    
    Nli1 = 12 + np.log10(interpolate.interp1d(t1['initial_mass'],t1['surface_li7'],\
            fill_value='extrapolate')(star_mass5))
   
    #convert teff and nli to EW   
    deltaEW = teff_nli_to_li(teff,Nli1) - teff_nli_to_li(teff,Nli5)
    deltaEW = np.clip(deltaEW,0,None)
    assert all(deltaEW >= 0), "EW must decrease monotonically"
    final_li = ngc2264_fit(const.BV) + deltaEW

    if (saveToFile):
        pickle.dump(final_li,open(join('data','mist_primordial_li.p'),'wb+'))
    return interpolate.interp1d(const.BV,final_li, fill_value='extrapolate')



def get_fit_BIC(bv_m,fits,dof):
    # Mamajek: Typical errors ... ∼0.1 dex (e.g. Henry et al. 1996; Paulson et al. 2002; White, Gabor, & Hillenbrand 2007)
    SIG = 0.1  
    num_stars = sum([len(bv_m[i][0]) for i in range(len(bv_m))])
    print("num stars", num_stars)
    
    log_L_hat = 0
    chi2_sum = 0
    
    for i in range(len(bv_m)):
        chi2 = prob.chi_sqr(bv_m[i][1],fits[i][0](bv_m[i][0]),SIG)
        chi2_sum += np.sum(chi2)
        log_L_hat += np.sum(-chi2/2)
    
    #print("log_L_hat= ", L_hat)
    print("chi2= ", chi2_sum)
    
    return np.log(num_stars)*dof - 2*log_L_hat
    

# refreshes everything saved to files
def main():
    bv_m,fits = readData.read_calcium(fromFile=False,saveToFile=True)
    _,res_arr = my_fits.get_fit_residuals(bv_m,fits,'calcium',None,li_range=None,
                linSpace=False,scale_by_std= False,vs_age_fit=True,zero_center=True)
    my_fits.fit_histogram('calcium',residual_arr=res_arr,fromFile=False,saveToFile=True)
    
    
        
    const = utils.init_constants('lithium')    
    bv_m, upper_lim, fits = readData.read_lithium(fromFile=False,saveToFile=True)
    
    my_fits.MIST_primordial_li(ngc2264_fit=fits[const.CLUSTER_NAMES.index('NGC2264')][0],fromFile=False, saveToFile=True)
    _,res_arr= my_fits.get_fit_residuals(bv_m,fits,'lithium',upper_lim,li_range=None,linSpace=False,
                                        vs_age_fit=True,zero_center=True)
    my_fits.fit_histogram('lithium',residual_arr=res_arr,fromFile=False,saveToFile=True)
    
    

if  __name__ == "__main__":
    main()










