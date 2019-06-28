import numpy as np
from scipy import integrate
import bisect
from scipy.stats import norm

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
