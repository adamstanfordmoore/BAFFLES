import numpy as np
import probability as prob
METAL_NAME = 'lithium'
CLUSTER_AGES = [5.5,24,43.7,85,130,200,240,600,700,4000]
CLUSTER_NAMES = ['NGC2264',r'$\beta$ Pic','IC2602',r'$\alpha$ Per','Pleiades',\
                 'M35','M34','Coma Ber','Hyades','M67']
NUM_STARS = [123, 37, 27, 60, 128, 82, 49, 13, 50, 40]
MARKERS = ['+','^','2','D','s','X','>','d','o','x']
COLORS = ['C0','C1','C8','C2','C3','C9','cornflowerblue','darkgreen','C4','C5']
PRIM_COLOR = 'lightcoral'
PRIM_MARKER = 'v'
BLDB_COLOR = 'darkmagenta'
BLDB_MARKER='8'
DOWN_ARROW = u'$\u2193$'

BV_RANGE = [.35,1.9] 
METAL_RANGE = [0.5,3.2]
METAL_RANGE_LIN = np.power(10,METAL_RANGE)

def inRange(bv,rhk):
    if not (BV_RANGE[0] <= bv <= BV_RANGE[1]):
        return False
    if not (METAL_RANGE[0] <= rhk <= METAL_RANGE[1]):
        return False
    return True

GALAXY_AGE = 13000 #Myr
BIN_SIZE = 10

BLDB_LOWER_LIM = 0.7 # B-V lower limit for using bldb point
FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

BV = np.linspace(BV_RANGE[0],BV_RANGE[1],1000) #the axis of the 2D arrays
BV_S = np.linspace(BV_RANGE[0],BV_RANGE[1],64) # spaced .025 apart 
AGE = np.logspace(0,np.log10(GALAXY_AGE),1000) #in units of Myr
METAL = prob.polyspace(0.5,10**METAL_RANGE[1],1000).reshape(1,1,1000) 

BV_UNCERTAINTY = 0.01
NUM_BV_POINTS = 15 #Number of points to represent measurement gaussian in baffles.likelihood
MEASURE_ERR = 15 #mA

#including new general piecewise li_vs_age 
DEFAULT_MEDIAN_GRID = "grids/median_li_061620.npy"

ZERO_LI = 0.5 #in log scale
PRIMORDIAL_LI_AGE = 1 #Myr
