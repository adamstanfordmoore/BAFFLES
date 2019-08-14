import numpy as np
import probability as prob
METAL_NAME = 'calcium'
CLUSTER_INDEX = [[0,9],[15,25],[9,15],[25,33],[33,46],[46,95],[181,198],[95,181],[198,274]]
CLUSTER_AGES = [10,16,24,45,85,130,500,625,4000]
CLUSTER_NAMES = ['Upper Scorpius','UCL+LCC',r'$\beta$ Pic','Tuc/Hor',r'$\alpha$ Per','Pleiades','UMa', 'Hyades','M67']

COLORS = ['C6','C7','C1','cornflowerblue','C2','C3','darkcyan','C4','C5']
MARKERS = ['H','*','^','4','D','s','P','o','x']

FMTS = [COLORS[i]+MARKERS[i] for i in range(len(MARKERS))]
#MARKERS = ['^','>','p','+','s','d','o','x']
#FMTS = ['C0s' for i in range(len(MARKERS))]
#FMTS = ['C0' + MARKERS[i] for i in range(len(MARKERS))]

#just main 4 clusters
#CLUSTER_AGES = [10,130,625,4000] #sco_cen,pleiades,hyades,m67
#CLUSTER_NAMES = ['Upper Sco','Pleiades','Hyades','M67']
#CLUSTER_INDEX = [[0,25],[46,95],[95,181],[198,274]]
#CLUSTER_INDEX = [[0,9,15,25],[46,95],[95,181],[198,274]] #no beta pic
#MARKERS = ['^','s','o','x']

#omitting pleiades
#CLUSTER_AGES = [10,16,24,85,500,625,4000]
#CLUSTER_NAMES = ['USco','UCL+LCC','Beta Pic', 'alpha Per','UMa', 'Hyades','M67']
#CLUSTER_INDEX = [[0,9],[15,25],[9,15],[33,46],[181,198],[95,181],[198,274]

#COLORS = ['C0','C4','C5','C6','C1','C7','C2','C3']
#COLORS = ['C0','C4','C2','C1','C5','C6','C7','C8']
#COLORS = ['C0','C1','C3','C6','C1','C7','C2','C3','C8']
#COLORS = ['C0','C1','C3','C6','C2','C7','C4','C5','C8']
#COLORS = ['C%s' % n for n in range(10)]

BV_RANGE = [.45,.9]
METAL_RANGE = [-5,-3.7]

def inRange(bv,rhk):
    if bv is not None and not (BV_RANGE[0] <= bv <= BV_RANGE[1]):
        return False
    if not (METAL_RANGE[0] <= rhk <= METAL_RANGE[1]):
        return False
    return True


GALAXY_AGE = 13000 #Myr 
BIN_SIZE = 10
DOWN_ARROW = u'$\u2193$'

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

BV = np.linspace(BV_RANGE[0],BV_RANGE[1],1000)#the axis of the 2D arrays
BV_S = np.array([.65]) #np.linspace(BV_RANGE[0],BV_RANGE[1],2)#lower sampling
METAL = prob.polyspace(METAL_RANGE[0],METAL_RANGE[1],1000) #the axis of the 2D arrays
AGE = np.logspace(0,np.log10(GALAXY_AGE),1000) #in units of Myr
#AGE = np.arange(1,GALAXY_AGE + 1,1)#np.logspace(0,4,1000) #in units of Myr
BV_UNCERTAINTY = .002
MEASURE_ERR = None #no default uncertainty in measurement 

#new numerical likelihood pdf from histogram fitting. still constant fits. polynomial fit at a given B-V like mamajek polynomial
DEFAULT_MEDIAN_GRID = "grids/median_rhk_081419.npy"
DEFAULT_SIGMA_GRID = "grids/sigma_rhk_081419.npy"

#new scatter from histogram fitting which equally weights each star.still constant fits. polynomial fit at a given B-V like mamajek polynomial
#DEFAULT_MEDIAN_GRID = "grids/median_rhk_030719.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_rhk_030719.npy"

#updated ages. constant fits/scatter. polynomial fit at a given B-V like mamajek polynomial
#DEFAULT_MEDIAN_GRID = "grids/median_rhk_103018.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_rhk_103018.npy"

#constant fits/scatter. polynomial fit at a given B-V like mamajek polynomial
#DEFAULT_MEDIAN_GRID = "grids/median_rhk_const.npy"
#DEFAULT_SIGMA_GRID = "grids/sigma_rhk_const.npy"
