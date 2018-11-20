import numpy as np
CLUSTER_INDEX = [[0,9],[15,25],[9,15],[33,46],[46,95],[181,198],[95,181],[198,274]]
CLUSTER_AGES = [10,16,24,85,130,500,625,4000]
CLUSTER_NAMES = ['USco','UCL+LCC','Beta Pic', 'alpha Per','Pleiades','UMa', 'Hyades','M67']

#just main 4 clusters
#CLUSTER_AGES = [10,130,625,4000] #sco_cen,pleiades,hyades,m67
#CLUSTER_NAMES = ['Upper Sco','Pleiades','Hyades','M67']
#CLUSTER_INDEX = [[0,25],[46,95],[95,181],[198,274]]
#CLUSTER_INDEX = [[0,9,15,25],[46,95],[95,181],[198,274]] #no beta pic

#omitting pleiades
#CLUSTER_AGES = [10,16,24,85,500,625,4000]
#CLUSTER_NAMES = ['USco','UCL+LCC','Beta Pic', 'alpha Per','UMa', 'Hyades','M67']
#CLUSTER_INDEX = [[0,9],[15,25],[9,15],[33,46],[181,198],[95,181],[198,274]

MARKERS = ['^','s','o','x']
#MARKERS = ['^','>','p','+','s','d','o','x']
#COLORS = ['C0','C4','C5','C6','C1','C7','C2','C3']
COLORS = ['C0','C4','C2','C1','C5','C6','C7','C8']

BV_RANGE = [.45,.9]
METAL_RANGE = [-5,-3.7]

GALAXY_AGE = 13001 #Myr with 1 for range purposes
BIN_SIZE = 10
DOWN_ARROW = u'$\u2193$'

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

BV = np.linspace(.45,.9,1000)#the axis of the 2D arrays
METAL = np.linspace(METAL_RANGE[0],METAL_RANGE[1],1000) #the axis of the 2D arrays
AGE = np.arange(1,GALAXY_AGE,1)#np.logspace(0,4,1000) #in units of Myr
BV_UNCERTAINTY = .002

#updated ages. constant fits/scatter. polynomial fit at a given B-V like mamajek polynomial
DEFAULT_MEDIAN_GRID = "grids/median_rhk_103018.npy"
DEFAULT_SIGMA_GRID = "grids/sigma_rhk_103018.npy"

#constant fits/scatter. polynomial fit at a given B-V like mamajek polynomial
#DEFAULT_MEDIAN_GRID = "grids/median_rhk_const.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_rhk_const.npy"
