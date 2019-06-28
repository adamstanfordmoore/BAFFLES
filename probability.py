import numpy as np
from scipy import integrate
import bisect
from scipy.stats import norm

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725] #[-2 sig,-1,mu,+1,+2]
UL_PROBS = [ 0.31731051,  0.04550026,  0.0026998,1] #[1-.68,1-.95,1-.99, maxAge]

#squared residuals divided by std**2 if given
def chi_sqr(x,mu,sig=1,total=False):
    chi2 = (x - mu)**2 / (sig**2)
    if (total):
        return np.sum(chi2)
    return chi2

def log_gaussian(x,mu,sig):
    chi2 = chi_sqr(x,mu,sig)
    return -np.log(sig*np.sqrt(2*np.pi)) - chi2/2

def lognorm(x,s):
    return 1/(s*x*np.sqrt(2*np.pi))*np.exp(-np.log(x)**2/(2*s**2))

#sig and mu can be numpy arrays
# x is either numpy array of same length or scalar
def gaussian(x,mu,sig):
    chi2 = chi_sqr(x,mu,sig)
    return 1/(sig*np.sqrt(2*np.pi)) * np.exp(-chi2/2)

def gaussian_cdf(x,mu,sig):
    return norm.cdf(x,mu,sig)

def polyspace(start,stop,num,power=2):
    x = np.arange(num)
    c = start
    a = (stop - start)/(num-1)**power
    return a*x**power + c

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

# Prior included for generality but is uniform for now
def prior(age):
    return 1
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
