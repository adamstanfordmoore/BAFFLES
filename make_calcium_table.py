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

def printStats(stats):
    arr = []
    for i,stat in enumerate(stats):
        FMT = "\\textbf{%d}"  if i==2 else "%d"
        arr.append(FMT % utils.round_sigs(stat,3))
    return arr

def main():
    ca_const = utils.init_constants('calcium')
    li_const = utils.init_constants('lithium')

    t = np.genfromtxt("data/nearbyStars_Boro_Saikia_2018.txt",delimiter=';',dtype=str,skip_header=57)

    baf_li = baffles.age_estimator('lithium')
    baf_ca = baffles.age_estimator('calcium')
    empty = ''

    f = open("nearbyStars_table",'w+')

    i = 0
    row = None
    count = 0

    while i < len(t):
        #print("I",i)
        if not utils.isFloat(t[i][10]) or not utils.isFloat(t[i][5]):
            i += 1
            continue
        elif row is None:
            row=copy.deepcopy(t[i])
            count += 1
            i += 1
            continue
        elif row[0] == t[i][0]:
            row[10] = float(t[i][10]) + float(row[10])
            row[5] = float(t[i][5]) + float(row[5])
            count += 1
            i += 1
            continue
        else:
            #print("HERE",row,t[i],count)
            row[10] = float(row[10])/count
            row[5] = float(row[5])/count
            count = 0

        #print(row)
        arr = []
        arr.append(row[3].strip())
        arr.append(row[0].strip())
        arr.append(row[1].strip())
        arr.append(row[4].strip())

        #if not utils.isFloat(row[10]) or not utils.isFloat(row[5]): continue
        bv = float(row[5])
        rhk = float(row[10])
        if not ca_const.inRange(bv,rhk):
            row = None
            continue
        arr.append("%.2f" % bv)
        arr.append("$%.3f$" % rhk)

        p_ca = baf_ca.get_posterior(bv,rhk,showPlot=False)
        arr += printStats(p_ca.stats)
        
        f.write(' & '.join(arr) + " \\\\")
        f.write('\n')

        row = None

    f.close()
    
if __name__ == '__main__':
    main()
