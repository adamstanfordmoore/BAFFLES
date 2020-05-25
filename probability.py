import numpy as np
from scipy import integrate
import bisect
from scipy.stats import norm
import matplotlib.pyplot as plt

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725] #[-2 sig,-1,mu,+1,+2]
UL_PROBS = [0.0026998,0.04550026,0.31731051,1] #[1-.99,1-.95,1-.68, maxAge]

#squared residuals divided by std**2 if given
def chi_sqr(x,mu,sig=1,total=False):
    if sig is None: sig = 1
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

# Space points according to the area of a guassian
# Divide gaussian into num equal area sections and put a point in the middle of each
# The result is denser sampling near the mean and less near the tails 
def gaussian_cdf_space(mu,sig,num,sig_lim=5):
    import fitting as my_fits
    x = np.linspace(mu - sig_lim*sig,mu + sig_lim*sig,300)
    cdf = gaussian_cdf(x,mu,sig)
    fit = my_fits.piecewise(cdf,x)
    temp = np.linspace(0,1,num+1) #splitting into equal area sections
    y = (temp[1:] + temp[:-1])/2 #making point the center of each section
    arr = fit(y)
    #plt.plot(x,gaussian(x,mu,sig))
    #plt.scatter(arr,gaussian(arr,mu,sig))
    #plt.show()
    return arr

# takes in densely sampled x,y and returns num sampled x
def desample(x,y,num):
    import fitting as my_fits
    ddy = np.abs(np.gradient(np.gradient(y)))
    f_arr = cdf(x,ddy)
    f_arr = f_arr*(x[-1] - x[0]) + x[0]
    f = my_fits.piecewise(f_arr,x)
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
def hist_cdf(vals):
    a,b = min(vals),max(vals)
    x = np.linspace(a - .05*(b-a),b + .05*(b-a),500)
    cdf = np.array([(vals < n).sum() for n in x],dtype='float')
    cdf /= cdf[-1]
    return x,cdf

#finds the x value with the largest y value
def mode(x,y):
    return x[np.argmax(y)]

#finds median age,ranges for 1,2 sigma as [-2 sigma,-1 sigma, median,+1 sigma,+2 sigma]
#finds ages from cdf values
def stats(age,y,upperLim=False): 
    import fitting as my_fits
    c = cdf(age,y)
    probs = UL_PROBS if upperLim else GAUSS_PROBS
    fit = my_fits.piecewise(c,age)
    return fit(probs)

# takes in a PDF given by age,y and a given age to compare to
# returns the percentile X such that given Age is within X %
def get_percentile(age,y,givenAge):
    import fitting as my_fits
    c = cdf(age,y)
    fit = my_fits.piecewise(age,c)
    return np.abs(fit(givenAge) - .5)*2*100

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
