import numpy as np
from scipy import integrate
import bisect
from scipy.stats import norm
from scipy.stats.mstats import gmean
import sys
import datetime

#returns true if sig or any element of sig is <= MIN_SIG
def negative_sig(sig):
    MIN_SIG = 1e-2
    return np.any(np.array(sig) < MIN_SIG)

def isFloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def hasNan(x):
    return np.any(np.isnan(x))

# takes in the logR'HK value and returns the age in units Myr
def getMamaAge(r):
    return np.power(10,-38.053 - 17.912*r - 1.6675*r*r)/1e6

def getMamaProductAge(r_arr):
    return gmean(getMamaAge(r_arr))

# Takes in age in units Myr and converts to logR'HK
def getMamaRHK(age):
    log_age = np.log10(np.array(age)*1e6)
    #Inverting Equation 3
    #a,b,c = -1.6675,-17.912,-38.053 - log_age
    #return (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
    #Equation 4
    return 8.94 - 4.849*log_age + .624*log_age**2 - .028 * log_age**3

def init_constants(metal):
    if (metal[0].lower() == 'c'):
        import ca_constants as const
    elif (metal[0].lower() == 'l'):
        import li_constants as const
    else:
        raise RuntimeError("No metal specified. Please enter lithium or calcium")
    return const

def progress_bar(frac,secondsLeft=None):
    sys.stdout.write('\r')
    sys.stdout.write("[%-25s] %d%% ETA: " % ('='*int(frac*100/4 - 1) + '>', frac*100) + \
            str(datetime.timedelta(seconds=int(secondsLeft))))
    if frac == 1:
        sys.stdout.write('\n')
    sys.stdout.flush()

def round_sigs(x,sigs=3):
    i = int(np.log10(x)+1)*-1 + sigs
    return np.round(x,i)

def float_sptype(sp):
    if sp[-1] == 'V' or sp[-1]== 'v':
        sp = sp[:-1]
    if '+' in sp:
        sp = sp.split('+')[0]
    S = 'OBAFGKMLTY'
    return S.index(sp[0])*10 + float(sp[1:])



