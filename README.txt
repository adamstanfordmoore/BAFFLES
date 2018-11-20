Adam Stanford-Moore,Eric Nielsen, Bruce Macintosh, Rob De Rosa
Stanford University Physics Department
8/28/18

This package implements BAFFLES: Bayesian Ages for Field LowEr-mass Stars, using measurements of calcium emission strength (logR'HK) and lithium abundance (log(EW/mA)).  This file contains a short decription of the process and header files for use.  

---------

The idea is to use Bayesian statistics to find a posterior age distribution for a star given its B-V and R’HK parameters (or logEW for lithium).
Mathematically we want P(age | R’HK,B-V), which using Bayes Equation becomes P(R’HK | age, B-V)*P(age).

(B-V)o: Star color corrected for dust
logR’HK: calcium emission strength related to relative flux from calcium emmision in chromosphere of stars - younger stars are spinning faster, stronger magnetic field, hotter chromosphere/corona, greater calcium emission
log(EW/mA): Lithium abundance measured by equivalent width of Li I 6708 A line, corrected for Fe I 6707.441. 

Illustrating method with Ca:
For each cluster I found a model for R’HK as a function of B-V (simple linear fit).
For each cluster I also found a model for scatter of R’HK from fit as a function of B-V. Most recently this was the constant standard deviation. 
For a given B-V value I can find the model R’HK and sigma_R’HK for each cluster (from my fit lines), and therefore can plot R’HK and sigma_R’HK against age (where each point is from a given cluster).

Therefore I made two large 2D grids of R’HK and sigma_R’HK, with axes of B-V and Age.  For each B-V value on the y axis I used my R’HK (or sigma_R’HK) vs age graphs from above.  Intermediate grids complete the first half the battle which is done beforehand.

Now for a given star with parameters (rhk,bv) we evaluate P(age | rhk,bv) = P(rhk | age, bv)*P(age) = P(rhk | age, bv) in linear space. so for an array of possible ages we compute the likelihood at each possible age: 1/sqrt(2pi)/sigma(age,bv) * e^(-chi^2/2) where chi^2 = (rhk - R’HK(age,bv))^2/(sigma(age,bv).  Where the functions f(age,bv) are the grids computed above.  To handle uncertainty in bv say given bv_o ± s, we do the mentioned computation weighted by the probability of getting a particular bv: 1/sqrt(2pi)s * e^((bv - bv_o)^2/2s).  P(age) is a uniform prior representign uniform star formation rate between 0 and 13 Gyrs ago.  

-------------

Files included in BAFFLES Package:
README.txt : this file
hd984.py : simple example utilizing baffles to make posteriors
baffles.py : core of the BAFFLES implementation.  includes age_estimator class and 
probability.py : functions for computing posterior confidence interval, gaussians, normalizing, and chi-squared
fitting.py :  functions for fitting to logR'HK vs B-V data (likewise logEW vs b-v).  maximizing log likelihood
plotting.py :  functions to plot posteriors, posterior products, and intermediate graphs
ca_constants.py : constants relating to calcium
li_constants.py : constants relating to lithium
data : directory containing the raw calcium and lithium data plus pre-made grids and lithium pickle files
readData.py : script for reading in lithium data and calcium data

______________

Baffles is written in python 2.7 with matplotlib version 2.2.2

baffles.py headers for user convenience.

# includes stats and array.
# stats -> 5 element array with age values for [.023,.159,.500,.841,.978] percent of cumulative posterior.  stats[2] = median.
# stats[3], stats[1] is 68% confidence interval.  stats[4],stats[0] is 95% confidence interval range
# array is the posterior array with one element for each age
class posterior:


# defines BAFFLES with a single metal either caclium or lithium
class age_estimator:

	#takes in a metal idicator either 'calcium' or 'lithium' denoting which method to use
        # option to input the grid_median and grid_sigma as arrays or as strings referencing saved .npy files
        def __init__(self,metal,grid_median=None,grid_sigma=None,default_grids=True):

	# option to input the grid_median and grid_sigma as arrays or as strings referencing saved .npy files
        def set_grids(self,grid_median,grid_sigma):

	"""
        Takes in bv the (B-V)o corrected color and the metallicity to return a posterior object.
        Metallicity: log(R'HK) if refering to calcium. log equivalent width per mA if lithium.
        pdfPage : optional object that if given will save pdf of posterior to it
	showPlot : boolean indicating whether to show the plot before saving it
	givenAge : number to be used to draw verticle red line on posterior plot with label 'Given Age'
	bv_uncertainty : uncertianity in B-V with default .002
	mamajekAge : optional number to be draw on plot with label 'Mamajek Age'
	Returns a posterior object with stats and array
	"""
        def get_posterior(self,bv,metallicity,pdfPage=None,showPlot=False,givenAge=None,bv_uncertainty=BV_UNCERTAINTY,mamajekAge=None):

	"""
	Plots a posterior product from two arrays of B-V and metallicity where each index represents a single star.
	bv_arr : array of individual stars' B-V values
	metallicity_arr : each star's logR'HK or logEW depending on with metal in use.  
	showStars : boolean whether to plot individual stars' posteriors
	title : optional title
	bv_errs : array with the uncertainties in B-V.  must be same length as bv_arr
	Returns posterior object of posterior product
	def posterior_product(self,bv_arr,metallicity_arr,pdfPage=None,showPlot=False,showStars=False,title=None,givenAge=None,bv_errs=None):

	#calculates and returns a 2D array of sigma b-v and age
	# optional savefile names
        def make_grids(self,fits,medianSavefile=None,sigmaSavefile=None,setAsDefaults=False):
	
	# returns the arrays holding the grids as grid_median, grid_sigma
	def get_grids(self):




___________________


plotting.py header file

# Plots the posterior y against age.
# starArray  is an array of posteriors for individual stars if plotting multiple posteriors on same graph in background. i.e. posterior product
def posterior(age,y,stat,title='Posterior Plot',pp=None,showPlot=False,starArray = [],givenAge=None,mamajekAge=None):

# plots the fits from fits onto bv_li data from each cluster
# bv_li : array of [B-v array,R'HK array] elements for each cluster
# fits : array of [mean fit function, scatter function] elements for each cluster. functions of B-V
# upper_lim : array of array for each cluster of boolean values representign if matellicity measurement is an upper limit
# e.g. [[F,F,F,F],[T,T,F],[T,F,F,F,T]]
# shadeScatter : boolean whether or not to shade the scatter around fit
# titles : optional array of titles for each cluster plot
def metal_vs_bv(bv_li,fits,metal,pdfPage=None,showPlots=False,upper_lim=None,titles=None):

# plots metal vs age at a given bv value
# fits : [[median func, scatter func] for each cluster]
def metal_vs_age(fits,metal,bv =.65,pdfPage=None,showPlots=False,title=None,shadeScatter=False):


_____________

#reads in data and generates fits.  returns bv_rhk, fits
# fromFile = True reads from pickle file
def read_calcium(fromFile=True,saveToFile=False):

#return bv_li, upper_lim, fits
def read_lithium(fromFile=True,saveToFile=False):

# returns the fits
def get_li_fits(bv_li,upper_lim_all):

_____________________________
