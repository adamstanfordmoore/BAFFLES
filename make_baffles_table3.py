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
import warnings



def abdor(annotate=False):
    import li_constants as const
    abdor_c = []
    abdor_l = []
    abdor_bverr = []
    abdor_lerr = []
    B,V,Name,SPT,REF = [],[],[],[],[]
    t = ascii.read('data/ab_dor_updated_err_bv_bib.csv',delimiter=',')
    for i,line in enumerate(t[1:]):
        if line[14] != '': continue 
        if not utils.isFloat(line[13]): 
            t[i+1,14] = 'OOR'
            continue
        bv = float(line[13])
        ew = float(line[3])
        if in_bounds(bv,ew,const):
            abdor_c.append(bv)
            abdor_l.append(ew)
            abdor_lerr.append(float(line[4]))
            B.append(line[10])
            V.append(line[11])
            SPT.append(line[12])
            REF.append(line[15])
            Name.append(line[0])
        else:
            t[i+1,14] = 'OOR'
            
    if annotate:
        np.savetxt("data/ab_dor_updated_err_bv_annotated.csv",t,delimiter=',',fmt='%s')
    
    
    abdor_c,abdor_l = np.array(abdor_c),np.log10(np.array(abdor_l))
    return abdor_c, abdor_l,abdor_lerr,B,V,SPT,REF,Name

def tuchor(annotate=False):
    import li_constants as const
    tuchor_c = []
    tuchor_l = []
    tuchor_bverr,tuchor_lerr = [],[]
    B,V,Name,SPT,REF = [],[],[],[],[]
    t = ascii.read('data/tuchor_updated_err_bv_bib.csv',delimiter=',')
    for i,line in enumerate(t[1:]):
        if line[14] != '' : continue 
        if not utils.isFloat(line[13]): 
            t[i+1,14] = 'OOR'
            continue
        bv = float(line[13])
        ew = float(line[3])
        if in_bounds(bv,ew,const):
            tuchor_c.append(bv)
            tuchor_l.append(ew)
            tuchor_lerr.append(float(line[4]))
            B.append(line[10])
            V.append(line[11])
            SPT.append(line[12])
            REF.append(line[15])
            Name.append(line[0])
        else:
            t[i+1,14] = 'OOR'
    
    if annotate:
        np.savetxt("data/tuchor_updated_err_bv_annotated.csv",t,delimiter=',',fmt='%s')
    
    tuchor_c,tuchor_l = np.array(tuchor_c),np.log10(np.array(tuchor_l))
    return tuchor_c,tuchor_l,tuchor_lerr,B,V,SPT,REF,Name
    
def schkolnik_betaPic(annotate=False):
    import li_constants as li_const
    t = np.genfromtxt("data/Shkolnik_2017_betaPic_bib.csv",delimiter=',',dtype=str)
    
    b,l,ul,names = [],[],[],[]
    B,V,SPT,REF,Name = [],[],[],[],[]
    for i,row in enumerate(t[1:]):
        if row[13] != '': continue
        if not utils.isFloat(row[1]) or not utils.isFloat(row[5]) \
                or row[12] != 'Y' or float(row[5]) <= 0: 
            t[i+1,13] = 'OOR'
            continue
        bv,ew = float(row[1]),np.log10(float(row[5]))
        if not li_const.inRange(bv,ew): 
            t[i+1,13] = 'OOR'
            continue
        b.append(bv)
        l.append(ew)
        ul.append(row[4]=='<')
        names.append(row[0])
        B.append(-999)
        V.append(-999)
        SPT.append(row[2])
        REF.append(row[16])
        Name.append(row[15])
    
    if annotate:
        np.savetxt("data/Shkolnik_2017_betaPic_annotated.csv",t,delimiter=',',fmt='%s')
    return b,l,ul,names,B,V,SPT,REF,Name

def mentuch2008_betaPic(annotate=False):
    import li_constants as li_const
    bp_c = []
    bp_l = []
    lim_bp = []
    names = []
    bp_err = []
    B,V,SPT,REF,Name = [],[],[],[],[]
    t = np.genfromtxt('data/beta_pic_updated_err_2MASS_bib.txt', delimiter='\t',dtype=str,skip_header=2)
    for i,line in enumerate(t):
        if line[11] != '': continue
        if in_bounds(float(line[4]),float(line[5]),li_const):
            bp_c.append(float(line[4]))
            bp_l.append(float(line[5]))
            lim_bp.append(False)
            names.append(line[10])
            bp_err.append(float(line[6]))
            B.append(line[2])
            V.append(line[3])
            SPT.append(line[9])
            REF.append(line[12])
            Name.append(line[0])
        else:
            t[i,11] = 'OOR'
    bp_c,bp_l = np.array(bp_c),np.log10(np.array(bp_l))
    
    if annotate:
        np.savetxt("data/beta_pic_updated_err_2MASS_annotated.txt",t,delimiter='\t',fmt='%s')
    return bp_c,bp_l,lim_bp,bp_err,names,B,V,SPT,REF,Name

def merged_betaPic():
    bp_c,bp_l,lim_bp,bp_err,names,B,V,SPT,REF,Name = mentuch2008_betaPic()
    bp_c2,bp_l2,lim_bp2,names2,B2,V2,SPT2,REF2,Name2 = schkolnik_betaPic()
    combined_err = [None]*len(names2)
    for i in range(len(names)):
        if names[i] in names2:
            j = names2.index(names[i])
            #print(names[i],bp_c[i],bp_c2[j],bp_l[i],bp_l2[j],lim_bp[i],lim_bp2[j])
            bv = (bp_c[i]+bp_c2[j])/2  #Not really averaging since b-v is the same for both
            bp_c2[j] = bv
            ew = None
            if lim_bp[i] and lim_bp2[j]:
                ew = min(bp_l[i],bp_l2[j])
                lim_bp2[j] = True
                #print("1")
            elif not lim_bp[i] and not lim_bp2[j]:
                
                ew = (bp_l[i] + bp_l2[j])/2
                lim_bp2[j] = False
                #print("2")
            else:
                ew = bp_l[i] if lim_bp2[j] else bp_l[i]
                lim_bp2[j] = False
                #print("3")
            bp_l2[j] = ew
            continue
        else:
            bp_c2.append(bp_c[i])
            bp_l2.append(bp_l[i])
            lim_bp2.append(lim_bp[i])
            names2.append(names[i])
            combined_err.append(bp_err[i])
            Name2.append(Name[i])
            REF2.append(REF[i])
            SPT2.append(SPT[i])
            V2.append(V[i])
            B2.append(B[i])
    return bp_c2,bp_l2,lim_bp2,combined_err,names2,B2,V2,SPT2,REF2,Name2    

#simple logic to check if (bv,l) is in bounds
#log represents if l is log(EW)
def in_bounds(bv,l,const,log=False):
    if (bv > const.BV_RANGE[0] and bv < const.BV_RANGE[1]):
        if (not log and (l > np.power(10,const.METAL_RANGE[0]) and l < np.power(10,const.METAL_RANGE[1]))):
            return True
        elif (log and (l > const.METAL_RANGE[0] and l < const.METAL_RANGE[1])):
            return True
    return False


def handleRefs(REF):
    star_bibs = [REF[i].split(';') for i in range(len(REF))]
    refs = np.unique(np.hstack(star_bibs)).tolist()
    
    bib_arr = []
    for i,bib in enumerate(refs):
        bib_arr.append('(%d) \citet{%s}' %(i+1,bib))
    latex_refs = ', '.join(bib_arr)
    
    numbers = []
    for row in star_bibs:
        inds = [str(refs.index(code) + 1) for code in row]
        numbers.append(','.join(inds))
    
    return numbers, latex_refs
    
def main():
    
    #c,l,err,B,V,SPT,REF,Name = [],[],[],[],[],[],[],[]
    ab_c, ab_l,ab_err,ab_B,ab_V,ab_SPT,ab_REF,ab_Name = abdor()
    t_c,t_l,t_lerr,t_B,t_V,t_SPT,t_REF,t_Name = tuchor()    
    
    #will need to keep track of references 
    #will need to add B,V for Schkolnik
    #will need 
    bp_c2,bp_l2,lim_bp2,combined_err,names2,B2,V2,SPT2,REF2,Name2  = merged_betaPic()

    c = np.concatenate((ab_c,t_c,bp_c2))
    B = np.concatenate((ab_B,t_B,B2))
    V = np.concatenate((ab_V,t_V,V2))
    SPT = np.concatenate((ab_SPT,t_SPT,SPT2))
    REF = np.concatenate((ab_REF,t_REF,REF2))
    Name = np.concatenate((ab_Name,t_Name,Name2))
    GROUP = np.concatenate((["AB Dor"]*len(ab_c),["Tuc/Hor"]*len(t_c),["$\\beta$ Pic"]*len(bp_c2)))
    
    Name = [x.strip().replace('V* ','') for x in Name]
    SPT = [x.strip() for x in SPT]


    ref_numbers,latex_refs = handleRefs(REF)
        
    headers = ["Name","SpT","Moving Group","B-V","ref"]
    
    
    f = open("baffles_table3_latex.txt",'w+')
    
    f.write(' & '.join(headers) + " \\\\")
    f.write('\n')
    f.write('\hline'+'\n')
    
    for i in range(len(c)):
        arr = [Name[i],SPT[i],GROUP[i],str(c[i]),ref_numbers[i]]
        f.write(' & '.join(arr) + " \\\\")
        f.write('\n')
    
    f.write('# ' + latex_refs)
    f.close()
    


if __name__ == '__main__':
    main()

