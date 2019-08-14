#!/usr/bin/env python2.7

from astropy.table import Table
from astropy.io import ascii
import numpy as np
np.seterr(divide = 'ignore')
import pickle
import fitting as my_fits
import utils
import copy
import baffles
import probability as prob

def printStats(stats):
    arr = []
    for i,stat in enumerate(stats):
        FMT = "\\textbf{%d}"  if i==2 else "%d"
        arr.append(FMT % utils.round_sigs(stat,3))
    return arr

def printHeader(title='BAFFLES Ages for stars in Nielsen 2010'):
    n_columns = 21
    #columns = ["Identifier","B-V","Sp. Type","logR'HK","LiEW (m\\AA)"]
    #columns += ["2.5% Age","16% Age","50% Age","84% Age","97.5% Age"]
    percents = ["2.5\\%","16\\%","\\textbf{50\\%}","84\\%","97.5\\%"]*3
    #columns += ["Group Membership"]


    print "\\begin{longrotatetable}"
    print "\\begin{deluxetable*}{" +  ''.join(['c' for _ in range(n_columns)]) + "}"
    print "\\tablewidth{0pt}"
    print "\\tabletypesize{\\tiny}"
    print "\\tablecaption{" + title + "}"
    print "\\startdata"
    print "Identifier & B-V & Sp. Type & logR'HK & LiEW & \\multicolumn{5}{c}{R'HK Age at Posterior CDF Value (Myr)} & \\multicolumn{5}{c}{LiEW Age at Posterior CDF Value (Myr)} & \\multicolumn{5}{c}{Final Age at Posterior CDF Value (Myr)} & Group Membership \\\\"
    print "& & & & (m\AA) &" + ' & '.join(percents) + " & \\\\"
    print "\hline"

def printFooter(label='nielsen_table'):
    print "\enddata"
    print "\label{table:" + label + "}"
    print "\end{deluxetable*}"
    print "\end{longrotatetable}"



def main():
    ca_const = utils.init_constants('calcium')
    li_const = utils.init_constants('lithium')

    t = np.genfromtxt('data/nielsen_2010_table2.csv',delimiter=',',dtype=str,skip_header=1)

    printHeader()

    baf_li = baffles.age_estimator('lithium')
    baf_ca = baffles.age_estimator('calcium')
    empty = ''
    
    bv_to_teff = my_fits.magic_table_convert('bv','teff')

    for row in t:
        arr = []
        arr.append(row[0])
        bv = None
        if utils.isFloat(row[1]):
            bv = float(row[1])
            arr.append("%.2f" % bv)
        else:
            arr.append(empty)
        
        arr.append(row[4].strip())
        
        p_ca,p_li = None,None
        if utils.isFloat(row[13]):
            rhk = float(row[13])
            arr.append(row[13])
            if ca_const.inRange(bv,rhk):
                p_ca = baf_ca.get_posterior(bv,rhk,showPlot=False)
        else:
            arr.append(empty)
       
        ew = None
        if utils.isFloat(row[7]):
            arr.append(row[7])
            ew = float(row[7])
        else:
            arr.append(empty)
        
        if bv is not None and ew is not None and li_const.inRange(bv,np.log10(ew)):
            p_li = baf_li.get_posterior(bv,np.log10(ew),showPlot=False)
        
        if p_ca is not None:
            arr += printStats(p_ca.stats)
        else:
            arr += [empty]*5

        if p_li is not None:
            arr += printStats(p_li.stats)
        else:
            arr += [empty]*5

        
        if p_ca is None and p_li is None:
            continue

        if p_ca is not None and p_li is not None:
            prod = p_ca.array * p_li.array
            prob.normalize(ca_const.AGE,prod)
            stats = prob.stats(ca_const.AGE,prod)
            arr += printStats(stats)
        else:
            arr += [empty]*5

        if row[16] != '':
            arr.append(row[16])
        else:
            arr.append(empty)
        
        print ' & '.join(arr) + "\\\\"

    printFooter() 
    
if __name__ == '__main__':
    main()
