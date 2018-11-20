import numpy as np
from scipy import integrate
import bisect

MIN_SIG = 1e-2
FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

#returns true if sig or any element of sig is <= MIN_SIG
def negative_sig(sig):
    return np.any(np.array(sig) < MIN_SIG)

def chi_sqr(mu,data,sig=1,total=False):
    chi2 = np.power(data - mu,2) / (sig*sig)
    if (total):
        return np.sum(chi2)
    return chi2

def log_gaussian(mu,sig,data):
    chi2 = chi_sqr(mu,data,sig)
    return -np.log(sig*np.sqrt(2*np.pi)) - chi2/2

def gaussian(mu,sig,data):
    chi2 = chi_sqr(mu,data,sig)
    return np.power(sig*np.sqrt(2*np.pi),-1) * np.exp(-chi2/2)

def normalize(x,y):
    area = np.trapz(y,x) 
    for i in range(len(y)):
        y[i] /= area

#scales y so that max(y) = height
def scale_to_height(y,height):
    scale = height/np.max(y)
    for i in range(len(y)):
        y[i] *= scale

#Cumulative Distribution Function 
def cdf(x,y): #depends on spacing 1 between age
    cum = integrate.cumtrapz(y,x) 
    return cum / cum[-1]

def stats(age,y): #finds median age,ranges for 1,2 sigma as [-2 sigma,-1 sigma, median,+1 sigma,+2 sigma]
    c = cdf(age,y)
    return [age[bisect.bisect_left(c,cum_prob)] for cum_prob in GAUSS_PROBS]

# Prior to be used in the case of using log-spaced ages
def log_prior(age):
    age_last = np.power(10,np.log10(age[-1]) + (np.log10(age[-1]) - np.log10(age[-2])))
    prior = [(age[i+1] - age[i])/(age[-1] - age[0]) for i in range(len(age)-1)]
    prior.append((age[i+1] - age[i])/(age[-1] - age[0]))
    return prior
