import numpy as np
import probability as prob
METAL_NAME = 'calcium'
CLUSTER_INDEX = [[0,9],[15,25],[9,15],[25,33],[33,46],[46,95],[181,198],[95,181],[198,274]]
CLUSTER_AGES = [10,16,24,45,85,130,500,700,4000]
CLUSTER_NAMES = ['Upper Sco','UCL+LCC',r'$\beta$ Pic','Tuc/Hor',r'$\alpha$ Per','Pleiades','UMa', 'Hyades','M67']
NUM_STARS = [8, 8, 6, 6, 12, 42, 10, 41, 70]


COLORS = ['C6','C7','C1','cornflowerblue','C2','C3','darkcyan','C4','C5']
MARKERS = ['H','*','^','4','D','s','P','o','x']

FMTS = [COLORS[i]+MARKERS[i] for i in range(len(MARKERS))]

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

FIVE_SIGMAS = 9.02e-07
GAUSS_PROBS = [.0227501,.158655,.5,.841345, .97725]

BV = np.linspace(BV_RANGE[0],BV_RANGE[1],1000)#the axis of the 2D arrays
BV_S = np.array([.65]) #lower sampling
METAL = prob.polyspace(METAL_RANGE[0],METAL_RANGE[1],1000) #the axis of the 2D arrays
AGE = np.logspace(0,np.log10(GALAXY_AGE),1000) #in units of Myr
BV_UNCERTAINTY = None #B-V not incorporated into fits
MEASURE_ERR = None #no default uncertainty in measurement 

#Fits as a function of age. 
#If path dividers are different on your operating system, run "python refresh.py"
DEFAULT_MEDIAN_GRID = "grids/median_rhk_061620.npy"

