#!/usr/bin/env python3

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
from make_nielsen_table import printStats,printHeader,printFooter


def main():
    ca_const = utils.init_constants('calcium')
    li_const = utils.init_constants('lithium')

    t = np.genfromtxt('data/brandt_2014_table.csv',delimiter=',',dtype=None,skip_header=2)
    
    printHeader("BAFFLES Ages for stars in Brandt 2014")

    f = open("brandt_2014_table",'w+')

    baf_li = baffles.age_estimator('lithium')
    baf_ca = baffles.age_estimator('calcium')
    empty = ''
    
    bv_to_teff = my_fits.magic_table_convert('bv','teff')

    for row in t:
        arr = []
        arr.append(row[0])
        bv = None
        if utils.isFloat(row[2]) and utils.isFloat(row[3]):
            bv = float(row[2]) - float(row[3])
            arr.append("%.2f" % bv)
        else:
            arr.append(empty)
        
        arr.append(row[4].strip().replace('_','\\_'))
        
        p_ca,p_li = None,None
        #rhk,li = row[7],row[9]
        if utils.isFloat(row[7]):
            rhk = float(row[7])
            arr.append('$'+row[7]+'$')
            if ca_const.inRange(bv,rhk):
                p_ca = baf_ca.get_posterior(bv,rhk,showPlot=False)
        else:
            arr.append(empty)
        
        ew = None
        if bv is not None and row[9].find('A') != -1:
            nli = float(row[9].strip()[-1])
            teff = bv_to_teff(bv)
            ew = 10** my_fits.teff_nli_to_li([teff],[nli])[0]
            arr.append("%d$^a$" % ew)
        elif utils.isFloat(row[9]):
            arr.append(row[9])
            ew = float(row[9])
        else:
            arr.append(empty)
        
        if bv is not None and ew is not None and li_const.inRange(bv,np.log10(ew)):
            p_li = baf_li.get_posterior(bv,ew,showPlot=False)
        
        
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


        if row[5] != '...':
            arr.append(row[5])
        else:
            arr.append(empty)

        f.write(' & '.join(arr) + " \\\\")
        f.write('\n')
    f.close()

    printFooter("brandt_table")

if __name__ == '__main__':
    main()
