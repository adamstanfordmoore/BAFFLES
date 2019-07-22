import numpy as np
from scipy import integrate
import bisect
from scipy.stats import norm
import fitting as my_fits
import matplotlib.pyplot as plt

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725] #[-2 sig,-1,mu,+1,+2]
UL_PROBS = [0.0026998,0.04550026,0.31731051,1] #[1-.99,1-.95,1-.68, maxAge]

#squared residuals divided by std**2 if given
def chi_sqr(x,mu,sig=1,total=False):
    if type(sig) == list:
        sig = np.array(sig) + 1e-10
    if type(mu) == list:
        mu = np.array(mu)
    if type(x) == list:
        x = np.array(x)

    chi2 = (x - mu)**2 / (sig**2)
    if (total):
        return np.sum(chi2)
    return chi2

def log_gaussian(x,mu,sig):
    chi2 = chi_sqr(x,mu,sig)
    return -np.log(sig*np.sqrt(2*np.pi)) - chi2/2

def lognorm(x,s):
    return 1/(s*x*np.sqrt(2*np.pi))*np.exp(-np.log10(x)**2/(2*s**2))

#sig and mu can be numpy arrays
# x is either numpy array of same length or scalar
def gaussian(x,mu=0,sig=1):
    chi2 = chi_sqr(x,mu,sig)
    return 1/(sig*np.sqrt(2*np.pi)) * np.exp(-chi2/2)

def gaussian_cdf(x,mu,sig):
    return norm.cdf(x,mu,sig)

def polyspace(start,stop,num,power=2):
    x = np.arange(num)
    c = start
    a = (stop - start)/(num-1)**power
    return a*x**power + c

# NOT FINISHED
# takes in densely sampled x,y and returns num sampled x
def desample(x,y,num):
    ddy = np.abs(np.gradient(np.gradient(y)))
    #ddy += 0.0001 ######### CHECK  #########
    f_arr = cdf(x,1/ddy)
    f_arr = f_arr*(x[-1] - x[0]) + x[0]
    f = my_fits.piecewise(x,f_arr)
    new_x = f(np.linspace(x[0],x[-1],num))
    new_y = my_fits.piecewise(x,y)(new_x)
    return new_x,new_y

#normalizes in-place
def normalize(x,y):
    area = np.trapz(y,x) 
    assert area > 0, "Invalid function to Normalize. Integral=%f" % area
    if area == 1: return y
    y[:] = y / area
    return y

#scales y so that max(y) = height
def scale_to_height(y,height):
    scale = height/np.max(y)
    y[:] = y*scale

#Cumulative Distribution Function 
def cdf(x,y): #depends on spacing 1 between age
    cum = integrate.cumtrapz(y,x,initial=0) 
    return cum / cum[-1]

#finds the x value with the largest y value
def mode(x,y):
    return x[np.argmax(y)]

#finds median age,ranges for 1,2 sigma as [-2 sigma,-1 sigma, median,+1 sigma,+2 sigma]
#finds ages from cdf
def stats(age,y,upperLim=False): 
    c = cdf(age,y)
    probs = UL_PROBS if upperLim else GAUSS_PROBS
    ages = []
    for cum_prob in probs:
        if cum_prob == 0: #handle edge case
            ages.append(age[0])
        else:
            i = bisect.bisect_left(c,cum_prob)
            a = (age[i] - age[i-1])/(c[i] - c[i-1]) * (cum_prob - c[i-1]) + age[i-1]
            ages.append(a)
    return ages

#calls func many times changing resample_args, keeping args constant
# returns product of calls
def resample(func,resample_args,args, sample_num=10,numIter=4):
    indices = np.arange(len(resample_args[0]))
    
    log_sum = 0
    for _ in range(numIter):
        inds = np.random.choice(indices,size=sample_num,replace=False)
        sampled_args = tuple(np.take(arr,inds) for arr in resample_args)
        argv = sampled_args + args
        log_sum += np.log(func(*argv))

    return np.exp(log_sum)





# Prior included for generality but is uniform for now
def prior(age,maxAge=None):
    agePrior = 1 if not maxAge else age <= maxAge
    return agePrior
    """
    if np.allclose(age[-1]-age[-2],age[1]-age[0]): #linear spacing
        return 1
    #otherwise non-uniform
    prior = (age[2:] - age[:-2])/2
    prior = np.insert(prior,0,age[1]-age[0])
    prior = np.append(prior,age[-1]-age[-2])
    normalize(age,prior)
    return prior
    """
