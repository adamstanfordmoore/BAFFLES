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
import sys

#If machine-readable, then drop $$ and add more sig figs
def printStats(stats, MR=False):
    arr = []
    for i,stat in enumerate(stats):
        FMT = "\\textbf{%d}"  if i==2 else "%d"
        if MR:
            arr.append('%g' % utils.round_sigs(stat,3))
        else:
            arr.append(FMT % utils.round_sigs(stat,3))
    return arr

def ra_dec(x):
    y = x.split()
    if len(y) == 6: #utils.isFloat(y[2]) and utils.isFloat(y[5]):
        y[3] = '$'+y[3] + '$'
        y[2] = '%.2f' % float(y[2])
        y[5] = '%.2f' % float(y[5])
    return ' '.join(y[:3]),' '.join(y[3:])

def average(list):
    sum,count = 0,0.0
    for x in list:
        if utils.isFloat(x):
            sum += float(x)
            count += 1
    return str(sum/count) if count > 0 else ''

# MR denotes machine_readable
def make_table(MR = False):
    ca_const = utils.init_constants('calcium')
    li_const = utils.init_constants('lithium')
    empty = ''
    """
    table = [] #[Object,RA,Dec,Sp Type,B-V,R'HK,Li EW,Source]
    #first read in all the 4 tables and create a single big table, which then I sort and merge

    t = np.genfromtxt('data/nielsen_2010_table2.csv',delimiter=',',dtype=str,skip_header=1)
    for row in t:
        if not utils.isFloat(row[1]) or not (.45 <= float(row[1]) <= 1.9): continue
        arr = []
        arr.append(row[21].strip())
        ra,dec = ra_dec(row[22])
        arr.append(ra)
        arr.append(dec)
        arr.append(row[4].strip())
        arr.append(row[1])
        arr.append(row[13])
        arr.append(row[7])
        arr.append("1")
        if arr[0] == '' or not (utils.isFloat(arr[5]) or utils.isFloat(arr[6])):
            continue
        table.append(arr)

    bv_to_teff = my_fits.magic_table_convert('bv','teff')
    t = np.genfromtxt('data/brandt_2014_table.csv',delimiter=',',dtype=str,skip_header=2)
    for row in t:
        bv = None
        if utils.isFloat(row[2]) and utils.isFloat(row[3]):
            bv = float(row[2]) - float(row[3])
        if bv is None or not (.45 <= bv <= 1.9): continue
        arr = []
        arr.append(row[14].strip())
        ra,dec = ra_dec(row[15])
        arr.append(ra)
        arr.append(dec)
        arr.append(row[4].strip())
        arr.append("%f" % bv)
        arr.append(row[7])
        if row[9].find('A') != -1:
            nli = float(row[9].split()[-1])
            teff = bv_to_teff(bv)
            ew = 10** my_fits.teff_nli_to_li([teff],[nli])[0]
            arr.append("%d" % ew)
        elif utils.isFloat(row[9]):
            arr.append(row[9])
        else:
            arr.append(empty)
        arr.append("2")
        if arr[0] == '' or not (utils.isFloat(arr[5]) or utils.isFloat(arr[6])):
            continue
        table.append(arr)


    t = np.genfromtxt("data/nearbyStars_Boro_Saikia_2018.txt",delimiter='\t',dtype=str,skip_header=58)
    for row in t:
        if not utils.isFloat(row[5]) or not (.45 <= float(row[5]) <= 1.9): continue
        arr = []
        arr.append(row[16].strip())
        ra,dec = ra_dec(row[17])
        arr.append(ra)
        arr.append(dec)
        arr.append(row[18].strip())
        arr.append(row[5])
        arr.append(row[10])
        arr.append(empty)
        arr.append("3")
        if arr[0] == '' or not (utils.isFloat(arr[5]) or utils.isFloat(arr[6])):
            continue
        table.append(arr)

    t = np.genfromtxt("data/guillot_2009_li_survey.txt",delimiter='\t',dtype=str,skip_header=77)
    for row in t:
        if not utils.isFloat(row[7]) or not (.45 <= float(row[7]) <= 1.9): continue
        arr = []
        arr.append(row[22].strip())
        ra,dec = ra_dec(row[23])
        arr.append(ra)
        arr.append(dec)
        arr.append(row[24].strip())
        arr.append(row[7])
        arr.append(empty)
        arr.append(row[16])
        arr.append("4")
        if arr[0] == '' or not (utils.isFloat(arr[5]) or utils.isFloat(arr[6])):
            continue
        table.append(arr)

    table = np.array(table)
    name_sorted = table[table[:,0].argsort()]

    thinned = [] #averaging b-v,measurements, sources as 1,4
    for name in set(name_sorted[:,0]):
        subset = name_sorted[name_sorted[:,0]==name]
        if len(subset) == 1:
            thinned.append(subset[0])
        else:
            arr = copy.deepcopy(subset[0])
            arr[4] = average(subset[:,4])
            arr[5] = average(subset[:,5])
            arr[6] = average(subset[:,6])
            x = list(set(subset[:,7]))
            x.sort()
            arr[7] = ','.join(x)
            thinned.append(arr)

    thinned = np.array(thinned)
    final_table = thinned[thinned[:,1].argsort()]
    np.save("final_table",final_table)
    exit()
    """

    final_table = np.load("data/merged_nielsen_brandt_saikia_guillot.npy")


    delimiterMR = ','

    baf_li = baffles.age_estimator('lithium')
    baf_ca = baffles.age_estimator('calcium')
    #[Object,RA,Dec,Sp Type,B-V,R'HK,Li EW,Source]
    f = open("baffles_table2_latex.txt",'w+') 
    fMR = open("baffles_table2.csv",'w+') 
    cdf = ['2.5%','16%','50%','84%','97.5%']
    column_head = ['Name','RA','Dec','Sp. Type','B-V',"logR'HK",'Li EW','Ref.']
    column_head += ["R'HK Age at CDF="+x for x in cdf]
    column_head += ["Li EW Age at CDF="+x for x in cdf]
    column_head += ["Final Age at CDF="+x for x in cdf]
    units = ['','h m s','h m s','','mags'," ",'mA','','','','','','','','','','','','','','','','']
    fMR.write(delimiterMR.join(column_head))
    fMR.write('\n')
    fMR.write(delimiterMR.join(units))
    fMR.write('\n')
    
    for row in final_table:
        arr = []
        arrMR = []
        arr += [x.replace('V* ','').replace('_','-') for x in row[0:4]] 
        arrMR += [x.replace('$','').replace('V* ','') for x in row[0:4]]
        
        bv = float(row[4])
        arr.append("%.2f" % bv) 
        arrMR.append("%.3g" % bv)        

        p_ca,p_li = None,None
        if utils.isFloat(row[5]):
            rhk = float(row[5])
            arr.append('$%.2f$' % rhk) 
            arrMR.append('%.3f' % rhk)
            if ca_const.inRange(bv,rhk):
                p_ca = baf_ca.get_posterior(bv,rhk,showPlot=False)
        else:
            arr.append(empty)
            arrMR.append(empty)
       
        ew = None
        if utils.isFloat(row[6]):
            ew = float(row[6])
            arr.append('%d' % ew)
            arrMR.append('%g' % ew)
        else:
            arr.append(empty)
            arrMR.append(empty)
        
        arr.append(row[7]) 
        arrMR.append(row[7].replace(',',';'))

        if bv is not None and ew is not None and ew > 0 and li_const.inRange(bv,np.log10(ew)):
            p_li = baf_li.get_posterior(bv,ew,showPlot=False)
        
        if p_ca is not None:
            arr += printStats(p_ca.stats)
            arrMR += printStats(p_ca.stats,MR=True)
        else:
            arr += [empty]*5
            arrMR += [empty]*5

        if p_li is not None:
            arr += printStats(p_li.stats)
            arrMR += printStats(p_li.stats,MR=True)
        else:
            arr += [empty]*5
            arrMR += [empty]*5

        
        if p_ca is None and p_li is None:
            continue

        if p_ca is not None and p_li is not None:
            prod = p_ca.array * p_li.array
            prob.normalize(ca_const.AGE,prod)
            stats = prob.stats(ca_const.AGE,prod)
            arr += printStats(stats)
            arrMR += printStats(stats,MR=True)
        elif p_ca is not None:
            arr += printStats(p_ca.stats)
            arrMR += printStats(p_ca.stats,MR=True)
        elif p_li is not None:
            arr += printStats(p_li.stats)
            arrMR += printStats(p_li.stats,MR=True)
        else:
            arr += [empty]*5
            arrMR += [empty]*5
             
        f.write(' & '.join(arr) + " \\\\")
        f.write('\n')
        fMR.write(delimiterMR.join(arrMR))   
        fMR.write('\n')
    f.close()
    fMR.close()


def main():
    #argv = sys.argv
    #make_table(MR = ('MR' in argv))
    make_table()




    
if __name__ == '__main__':
    main()
