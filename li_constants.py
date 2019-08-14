import numpy as np
import probability as prob
METAL_NAME = 'lithium'
#CLUSTER_AGES = [2,32,32,34.5,36,44,69,590]
#CLUSTER_NAMES = ['Cha I','IC4665','IC2602','NGC2451B','NGC2547','IC2391','NGC2516','NGC6633']
#CLUSTER_NAMES = ['Cha I','IC4665','IC2602','NGC2451A','NGC2451B','NGC2547','IC2391','NGC2516','NGC6633']
#CLUSTER_NAMES = ['Cha I','IC4665','IC2602','NGC2451B','NGC2547','IC2391','NGC2516']

CLUSTER_AGES = [5.5,24,43.7,85,130,200,240,600,625,4000]
CLUSTER_NAMES = ['NGC2264',r'$\beta$ Pic','IC2602',r'$\alpha$ Per','Pleiades',\
                 'M35','M34','Coma','Hyades','M67']

MARKERS = ['+','^','2','D','s','X','>','d','o','x']
COLORS = ['C0','C1','C8','C2','C3','C9','cornflowerblue','darkgreen','C4','C5']
PRIM_COLOR = 'lightcoral'
PRIM_MARKER = 'v'
BLDB_COLOR = 'darkmagenta'
BLDB_MARKER='8'
DOWN_ARROW = u'$\u2193$'


#MARKERS = ['^','>','p','+','s','d','o','x','<','P','*','D','H','o']
#MARKERS = ['o' for _ in range(len(CLUSTER_NAMES)+2)]
#COLORS = ['C%d' % i for i in np.arange(len(CLUSTER_NAMES)+2) % 10]
#COLORS = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9', '#1f77b4','#ff7f0e']
#COLORS = ['C0','C1','C2','C3','C4','C5','C6','C7','C4','C5']

BV_RANGE = [.35,1.9] #BV_RANGE = [.24,2.3]
#BV_RANGE = [.35,1.4] #BV_RANGE = [.24,2.3]
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
#AGE = np.arange(1,GALAXY_AGE + 1,1) #in units of Myr
METAL = prob.polyspace(0.5,10**METAL_RANGE[1],1000).reshape(1,1,1000) 
#METAL = np.logspace(-0.3,METAL_RANGE[1],1000).reshape(1,1,1000) 

BV_UNCERTAINTY = 0.01
NUM_BV_POINTS = 15 #Number of points to represent measurement gaussian in baffles.likelihood
MEASURE_ERR = 15 #mA


#including new general piecewise li_vs_age 
DEFAULT_MEDIAN_GRID = "grids/median_li_081419.npy"
DEFAULT_SIGMA_GRID = "grids/sigma_li_081419.npy"

#including prim_li, all clusters, new sigma scatter + measurement,\
  #constrained polynomial plus patch from bv .75 to .95
#DEFAULT_MEDIAN_GRID = "grids/median_li_030719.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_li_030719.npy"

#including BLDB. polnomial fits to li vs bv. interpolation at given B-V
#DEFAULT_MEDIAN_GRID = "grids/median_li_103018.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_li_103018.npy"

#polnomial fits to li vs bv. interpolation at given B-V
#DEFAULT_MEDIAN_GRID = "grids/median_li.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_li.npy"

ZERO_LI = 0.5 #in log scale
PRIMORDIAL_LI_AGE = 1 #Myr
