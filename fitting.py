"""
Module giving necessary fitting functions
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
import copy
import bisect
from scipy import integrate
from scipy.optimize import minimize
import probability as prob

BIN_SIZE = 10
MIN_PER_BIN = 4
MIN_LENGTH = .15 #minimum B-V length of piecewise segment
DOWN_ARROW = u'$\u2193$'
MEASUREMENT_LI_SCATTER = 10 #mA in linear space from measurement error

#limits vertex from moving to the right of lim
def constrained_poly_fit(x,y,lim=0):
    guess = np.polyfit(x,y,2)
    res = minimize(constrained_poly_minimizer,guess, args=(x,y,lim),method='Nelder-Mead')
    if (not res.success):
        print "Unsuccessful minimizing of constrained polynomial function, check initial guess"
        print res.x

    return np.poly1d(res.x)

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
    

def constrained_poly_minimizer(params,x,y,lim):
    #if max is greater than lim return a large number
    vert = -params[1]/(2*params[0])
    #vert = right_cubic_root(params)
    if (vert > lim):
        return 1000000*(vert - lim)
    y_model = np.poly1d(params)(x)
    return prob.chi_sqr(y_model,y,total=True)

#computes fit of age and lithium depletion boundary (well bv at which lithium goes to zero)
# returns function such that giving it a b-v value returns oldest age it could be.  
def ldb_fit(fits):
    bv_at_zero_li,ages,cluster_names = [],[],[]
    import li_constants as const
    for c in range(len(fits)):
        li = fits[c][0](const.BV)
        i = bisect.bisect_left(const.BV,.65) #only interested in zero crossing at bv > .65
        while (i < len(li) and li[i] > const.ZERO_LI):
            i += 1 #can change to 5 for slightly more speed
        if (i != len(li)):
            bv_at_zero_li.append(const.BV[i])
            ages.append(const.CLUSTER_AGES[c])
            cluster_names.append(const.CLUSTER_NAMES[c])
   
    """
    pp = PdfPages('ldb_vs_age.pdf')
    plt.title('Lithium depletion boundary over time')
    plt.xlabel('B-V',size=18)
    plt.ylabel('Age (Myr)',size=18)
    plt.yscale('log')
    """
    fit = poly_fit(bv_at_zero_li,np.log10(ages),1)
    def ldb_age(x):
        return np.power(10,fit(x))
    """ 
    for c in range(len(ages)):
        plt.scatter(bv_at_zero_li[c],ages[c],label=cluster_names[c])
    plt.plot(const.BV,ldb_age(const.BV))
    plt.legend()
    pp.savefig()
    plt.show()
    pp.close()
    """
    return ldb_age

def ldb_scatter(fits):
    bv_at_zero_li,ages,cluster_names = [],[],[]
    import li_constants as const
    for c in range(len(fits)):
        li = fits[c][0](const.BV)
        i = bisect.bisect_left(const.BV,.65) #only interested in zero crossing at bv > .65
        while (i < len(li) and li[i] > const.ZERO_LI):
            i += 1 #can change to 5 for slightly more speed
        if (i != len(li)):
            bv_at_zero_li.append(const.BV[i])
            ages.append(const.CLUSTER_AGES[c])
            cluster_names.append(const.CLUSTER_NAMES[c])
    fit = poly_fit(bv_at_zero_li,np.log10(ages),1)
    def ldb_age(x):
        return np.power(10,fit(x))
    return np.std(residuals(bv_at_zero_li,np.log10(ages),ldb_age))

#returns the 1d polynomial
def linear_fit(x,y,scatter=False):
    return poly_fit(x,y,1,scatter=scatter)

# finds and returns n dimensional polynomial fit.  If scatter=True than it returns [median fit,scatter fit]
def poly_fit(x,y,n=2,upper_lim=None,scatter=False):
    if (upper_lim):
        return minimize_polynomial(x,y,n,upper_lim)
    fit = np.poly1d(np.polyfit(x,y,n))
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


def poly_minimizer(params,c,r,upper_lim):
    sig = params[-1]
    if (prob.negative_sig(sig)):
        return np.inf
    r_model = np.poly1d(params[0:-1])(c)
    sig_model = log_scatter_with_linear_offset(sig,np.poly1d(params[0:-1]))(c)
    #sig_model = constant_fit(sig)[0](c)
       
    return inverse_log_likelihood(r_model,sig_model,r,upper_lim)

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
    upper_log = np.log10(linEW + sig_tot_lin)
    lower = linEW - sig_tot_lin
    lower[lower <=0] = 1

    
    sig_tot_log = (np.log10(linEW + sig_tot_lin) - np.log10(lower))/2
    return interpolate.interp1d(const.BV,sig_tot_log, fill_value='extrapolate')


#x_method: determines how x_coordinates are selected.
#   'even' means evenly spaced in the x dimension from min(x) to max(x)
#   'bin' means spaced evenly with number per bin at least 'bins' (the next parameter) in size
#   'free' means the x_coordinates are free parameters with only the min and max fixed by the range of x 
#returns piecewise function
def pwise_fit(x,y,segments=-1,upper_lim=None, guess_fit=None, x_method='even', bin_size=MIN_PER_BIN):
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
            x_coords = x_coords_smart_binning(x,bin_size)   
        else:
            x_coords = x_even(x,segments)
        params += x_coords[1:-1]
    elif (x_method[0] == 'b'):
        x_coords = x_coords_smart_binning(x,bin_size)
        
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
    def fit(x): #returns inner function
            return piecewise(x,x_val,y_val)
    def scatter_fit(x):
            return step(x,x_val[1:-1],sig)
    return [fit, scatter_fit] #,res.success

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

def compare_fits(params,argv): #c,r,upper_lim):
    x_coords,c,r,upper_lim = argv[0],argv[1],argv[2],argv[3]
    x,y,sig = get_x_y_sig(params,x_coords)
    r_model = piecewise(c,x,y)
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

def piecewise_minimizer(params,x_coords,c,r,upper_lim):
    x,y,sig = get_x_y_sig(params,x_coords)
    if (prob.negative_sig(sig)):
        return np.inf
        if (len(x) != len(y) or len(sig) != len(x) - 1):
                print params
                print x,y,sig
                raise RuntimeError("Error in inverse_log_likelihood")
        r_model = piecewise(c,x,y)
        sig_model = step(c,x[1:-1],sig)
        return inverse_log_likelihood(r_model,sig_model,r,upper_lim)

#minimizing the negative of the log_likelihood
def inverse_log_likelihood(r_model,sig_model,r,upper_lim):
    total_sum = 0
    for i in range(len(r)):
            if (upper_lim[i]):
                    r_range = np.linspace(r_model[i] - 4*sig_model[i],r[i],100)
                    total_sum += np.log(np.trapz(prob.gaussian(r_model[i],sig_model[i],r_range),r_range))
            else:
                    total_sum += prob.log_gaussian(r_model[i],sig_model[i],r[i])
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
def piecewise(x,x_locs,y_locs):
        return interpolate.interp1d(x_locs,y_locs, fill_value='extrapolate')(x)

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
    
#enforces bins to have no fewer than bi_size elements and
# have length no smaller than MIN_LENGTH
def x_coords_smart_binning(c,bin_size):
    a = copy.deepcopy(c)
    a.sort()
    x_coords = [a[0]]
    i = bin_size
    while (i < len(a)):
        new_x = a[i]#np.mean(a[i-1:i+1]) #average of two adjacent star's x_coords
        if (new_x - x_coords[-1] < MIN_LENGTH):
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

