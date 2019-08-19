# BAFFLES: Bayesian Ages for Lower-mass Field Stars

This package computes age posteriors for field stars from measurements of R'HK calcium
emission and/or B-V color and lithium equivalent width absorption.  For calcium emission 
our method is calibrated to stars with $B-V$ between 0.45 and 0.9 ($\sim$F6-K2) and 
log($R'_{HK}) between -3.7 and -5.  For lithium we have calibrated BAFFLES to stars with 
$B-V$ between 0.35 and 1.9 ($\sim$F2-M5) and EW Li between 3.2 and 1500 m\AA.  See 
the paper Stanford-Moore et al, 2020.      

### Requirements

Python 3
astropy (https://docs.astropy.org/en/stable/install.html#installing-astropy)

## Authors

* Adam Stanford-Moore
* Eric Nielsen
* Rob De Rosa
* Bruce Macintosh

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


## Command Line Usage Examples

Quick usage from the command line. Let's find an age for the sun using B-V of 0.65 and 
logR'HK = -4.908 (Mamajek & Hillenbrand 2008). By default it saves a pdf.
 
'./baffles.py -bmv 0.65 -rhk -4.906 -showPlot'

Now age of HR 2562 using B-V=.45 ± .02, log(R'HK) = -4.55 (Gray 2006), and lithium EW 
of 21 ± 5 (Mesa el al 2018). "-ul" would denote an upper limit.  "-s" will 
save a text file of the posterior. "-noPlot" will suppress plotting. -maxAge 10000 will 
constrain the prior on age to be uniform out to 10 Gyr. "-li\_err" allows input of uncertainty
on EW Li, and "-bv\_err" uncertainty on B-V. The following command will determine the age
using calcium and lithium separately and then find the combined posterior product.    

'./baffles.py -bmv 0.45 -bmv\_err .02 -rhk -4.55 -li 21 -li\_err 5 -showPlot'

Type "./baffles.py -help" to learn more.

To directly import baffles and use the module in a python script see "usage_example.py"


## Package Files

baffles.py    : main package that computes posteriors and final ages. makes and stores grids
of mean indicator as functions of age and B-V

fitting.py    : fitting functions for mean as a function of age and B-V
plotting.py   : plotting functions to display posteriors, data, and fits
probability.py   : assortment of statistical functions like finding gaussian PDF/CDF
readData.py    : reads in lithium and calcium data from data directory and returns data
                and indicator vs B-V fits
refresh.py   : updates any stored fits like mean activity as functions of age/B-V
utils.py    : extra helper functions
ca\_constants.py  :constants related to calcium
li\_constants.py  :constants related to lithium
usage\_example.py  : example script using the baffles package
paper\_plots\_li.py  : more advanced plotting examples for lithium used to create the plots 
in stanford_moore et al 2019.
paper\_plots\_ca.py   : more advanced plotting for calcium
data directory   :   contains files with calcium/lithium data and .p pickle files with
                    saved arrays and indicator vs B-V
grids directory   : contains saved grids of mean/sigma indicator values as functions of age/B-V







