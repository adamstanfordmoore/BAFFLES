import numpy as np
CLUSTER_AGES = [5.5,43.7,130,240,625,4000]
CLUSTER_NAMES = ['NGC2264','IC2602','Pleiades','M34','Hyades','M67']
MARKERS = ['^','>','p','+','s','d','o','x']

BV_RANGE = [.24,2.3]
METAL_RANGE = [0.5,3.2]
COLORS = ['C0','C1','C2','C4','C5','C6','C7','C8']
GALAXY_AGE = 13001 #Myr with 1 for range purposes
BIN_SIZE = 10
DOWN_ARROW = u'$\u2193$'

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

BV = np.linspace(BV_RANGE[0],BV_RANGE[1],1000) #the axis of the 2D arrays
AGE = np.arange(1,GALAXY_AGE,1)#np.logspace(0,4,1000) #in units of Myr
METAL = np.linspace(METAL_RANGE[0],METAL_RANGE[1],1000)

#including BLDB. polnomial fits to li vs bv. interpolation at given B-V
DEFAULT_MEDIAN_GRID = "grids/median_li_103018.npy"
DEFAULT_SIGMA_GRID = "grids/sigma_li_103018.npy"

#polnomial fits to li vs bv. interpolation at given B-V
#DEFAULT_MEDIAN_GRID = "grids/median_li.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_li.npy"

ZERO_LI = 0.5 #in log scale
