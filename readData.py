#!/usr/bin/env python2.7

"""
11/28/18
Imports the lithium data to three lists:
bv_li has an array for each cluster of [np.arra([B-V]),np.arrya([LogEW])]
bv_li = [[[cluster 1 B-V],[cluster1 LogEW]],[[cluster 1 B-V],[cluster1 LogEW]],....]
upper_lim has an array for each cluster with a boolean for each data point indicating if point is an upper limit
upper_lim [[T,F,T,F,F,F],[F,F,F,F,F]....]
fits stored as [median function,scatter function]
"""
import warnings
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import pickle
import fitting as my_fits

#reads in data and generates fits
def read_calcium(fromFile=True,saveToFile=False,fit_degree=0):
    if (fromFile):
        bv_rhk = pickle.load(open('data/bv_rhk.p','rb'))
        fits = pickle.load(open('data/ca_fits.p','rb'))
        return bv_rhk,fits
    import ca_constants as const
    warnings.simplefilter('ignore', UserWarning)
    t = Table.read('data/mamajek_table_5.fits')
    fits = []
    bv_rhk = []
    cluster_index = const.CLUSTER_INDEX
    for i in range(len(cluster_index)):
        c,r = [],[]
        if (len(cluster_index[i]) > 2):
            for j in range(0,len(cluster_index[i]),2):
                c += list(t["__B-V_0"][cluster_index[i][j]:cluster_index[i][j+1]])
                r += list(t["logR_HK"][cluster_index[i][j]:cluster_index[i][j+1]])
            c = np.array(c).astype(np.float)
            r = np.array(r).astype(np.float)
        else:
            c = np.array(t["__B-V_0"][cluster_index[i][0]:cluster_index[i][1]]).astype(np.float)
            r = np.array(t["logR_HK"][cluster_index[i][0]:cluster_index[i][1]]).astype(np.float)
        bv_rhk.append([c,r])
        #fits.append(my_fits.constant_fit(r))
        fits.append(my_fits.poly_fit(c,r,n=fit_degree,scatter=True))
    if (saveToFile):
        pickle.dump(bv_rhk,open('data/bv_rhk.p','wb'))
        pickle.dump(fits,open('data/ca_fits.p','wb'))
    return bv_rhk,fits

def get_li_fits(bv_li,upper_lim_all):
    import li_constants as const
    fits = []
    for i in range(len(bv_li)):
        poly_order = 2
        #if (const.CLUSTER_NAMES[i] in ['UMa','NGC3680']): #these cluster need linear fits
        #    poly_order = 1
        fit = my_fits.poly_fit(bv_li[i][0],bv_li[i][1],poly_order,upper_lim_all[i])
        #fit = my_fits.pwise_fit(bv_li[i][0],bv_li[i][1],upper_lim=upper_lim_all[i],guess_fit=my_fits.poly_fit(bv_li[i][0],bv_li[i][1]),x_method='bin')
        fits.append(fit)
    return fits

#simple logic to check if (bv,l) is in bounds
#log represents if l is log(EW)
def in_bounds(bv,l,const,log=False):
    if (bv > const.BV_RANGE[0] and bv < const.BV_RANGE[1]):
        if (not log and (l > np.power(10,const.METAL_RANGE[0]) and l < np.power(10,const.METAL_RANGE[1]))):
            return True
        elif (log and (l > const.METAL_RANGE[0] and l < const.METAL_RANGE[1])):
            return True
    return False

def read_lithium(fromFile=True,saveToFile=False):
    if (fromFile):
        bv_li = pickle.load(open('data/bv_li.p','rb'))
        upper_lim = pickle.load(open('data/upper_lim.p','rb'))
        fits = pickle.load(open('data/li_fits.p','rb'))
        return bv_li,upper_lim,fits

    import li_constants as const
    bv_li = []
    upper_lim = []
    t = ascii.read('data/ngc2264_lithium_bv.csv',delimiter=';')
    ngc2264_c = []
    ngc2264_l = []
    for line in t[2:]:  
        if (float(line[6]) != -999.0 and in_bounds(float(line[7]),1000*float(line[6]),const)):
            ngc2264_c.append(float(line[7]))
            ngc2264_l.append(1000*float(line[6]))

    t = ascii.read('data/ngc2264_lithium2.txt', delimiter=',')
    for line in t:
        if in_bounds(line[5],line[3],const):
            ngc2264_c.append(line[5])
            ngc2264_l.append(line[3]) 
    
    ngc2264_c,ngc2264_l = np.array(ngc2264_c),np.log10(np.array(ngc2264_l))
    ind = np.argmax(ngc2264_c)
    ngc2264_c = np.delete(ngc2264_c,ind)
    ngc2264_l = np.delete(ngc2264_l,ind)
    bv_li.append([ngc2264_c,ngc2264_l])
    upper_lim.append([False]*len(ngc2264_c))
    
    ic2602_c = []
    ic2602_l = []
    t = ascii.read('data/ic2602_lithium.txt', delimiter=',')
    for line in t:
            if in_bounds(line[1],line[3],const):
                    ic2602_c.append(line[1])
                    ic2602_l.append(line[3])
    ic2602_c,ic2602_l = np.array(ic2602_c),np.log10(np.array(ic2602_l))
    bv_li.append([ic2602_c,ic2602_l])
    upper_lim.append([False]*len(ic2602_c))
    
    pleiades_c = []
    pleiades_l = []
    lim_p = []
    t = ascii.read('data/pleiades_lithium.tsv', delimiter=';')
    for line in t[2:]:
        if in_bounds(float(line[3]),10*float(line[-7]),const):
                pleiades_c.append(float(line[3]))
                pleiades_l.append(10*float(line[-7]))
        if (line[-8] == '<'):
            lim_p.append(True)
        else:
            lim_p.append(False)
    pleiades_c,pleiades_l = np.array(pleiades_c),np.log10(np.array(pleiades_l))
    bv_li.append([pleiades_c,pleiades_l])
    upper_lim.append(lim_p)
    
    m34_c = []
    m34_l = []
    lim_m34 = []
    t = ascii.read('data/m34_lithium.txt', delimiter=',')
    for line in t:
        if (float(line[3]) < 0):
            if (in_bounds(float(line[1]),-float(line[3]),const)):
                m34_c.append(float(line[1]))
                m34_l.append(-float(line[3]))
                lim_m34.append(True)    
        elif (in_bounds(float(line[1]),float(line[3]),const)):
            m34_c.append(float(line[1]))
            m34_l.append(float(line[3]))
            lim_m34.append(False)
    m34_c,m34_l = np.array(m34_c),np.log10(np.array(m34_l))
    bv_li.append([m34_c,m34_l])
    upper_lim.append(lim_m34) 

    hyades_c = []
    hyades_l = []
    lim_h = []
    t = ascii.read('data/hyades_lithium.tsv', delimiter=';')
    for line in t[2:]:
        if in_bounds(float(line[5]),float(line[10]),const):
                hyades_c.append(float(line[5]))
                hyades_l.append(float(line[10]))
        if (line[9] == '<'):
                lim_h.append(True)
        else:
                lim_h.append(False)
    hyades_c,hyades_l = np.array(hyades_c),np.log10(np.array(hyades_l))
    bv_li.append([hyades_c,hyades_l])
    upper_lim.append(lim_h)

    m67_c = []
    m67_l = []
    lim_m67 = []
    t = ascii.read('data/m67_lithium.txt', delimiter=',')
    for line in t:
        if (float(line[4]) < 0):
            if in_bounds(float(line[2]),-float(line[4]),const):
                m67_c.append(float(line[2]))
                m67_l.append(-float(line[4])) 
                lim_m67.append(True)
        elif (in_bounds(line[2],line[4],const)):
            m67_c.append(line[2])
            m67_l.append(line[4])
            lim_m67.append(False)
        
    m67_c,m67_l = np.array(m67_c),np.log10(np.array(m67_l))
    bv_li.append([m67_c,m67_l])
    upper_lim.append(lim_m67)

    aper_c = []
    aper_l = []
    lim_aper = []
    t = ascii.read('data/alpha_per.txt', delimiter='\t')
    for line in t:
        if in_bounds(float(line[2]),float(line[4]),const):
            aper_c.append(float(line[2]))
            aper_l.append(float(line[4])) 
            lim_aper.append(False)
    aper_c,aper_l = np.array(aper_c),np.log10(np.array(aper_l))
    bv_li.append([aper_c,aper_l])
    upper_lim.append(lim_aper)

    coma_c = []
    coma_l = []
    lim_coma = []
    t = ascii.read('data/coma_berenices.txt', delimiter=',')
    for line in t:
        if in_bounds(float(line[2]),float(line[4]),const):
            coma_c.append(float(line[2]))
            coma_l.append(float(line[4]))
            lim_coma.append(line[-1]==1)
    coma_c,coma_l = np.array(coma_c),np.log10(np.array(coma_l))
    bv_li.append([coma_c,coma_l])
    upper_lim.append(lim_coma)

    m35_c = []
    m35_l = []
    lim_m35 = []
    t = ascii.read('data/M35_data.txt')
    for line in t:
        if in_bounds(float(line[1]),float(line[2]),const):
            m35_c.append(float(line[1]))
            m35_l.append(float(line[2]))
            lim_m35.append(line[-1]!=0)
    m35_c,m35_l = np.array(m35_c),np.log10(np.array(m35_l))
    bv_li.append([m35_c,m35_l])
    upper_lim.append(lim_m35)
    
    """
    uma_c = []
    uma_l = []
    lim_uma = []
    t = ascii.read('data/UMa.csv', delimiter=',')
    for line in t:
        if in_bounds(float(line[2]),float(line[3]),const):
            uma_c.append(float(line[2]))
            uma_l.append(float(line[3]))
            lim_uma.append(line[-1]==1)
    uma_c,uma_l = np.array(uma_c),np.log10(np.array(uma_l))
    bv_li.append([uma_c,uma_l])
    upper_lim.append(lim_uma)

    ngc3680_c = []
    ngc3680_l = []
    lim_ngc3680 = []
    t = ascii.read('data/ngc_3680.csv',delimiter=',')
    for line in t[1:]:
        Te = float(line[2])
        c = 8575 - Te
        BV = (5384.091 - np.sqrt(5384.091*5384.091 - 4*1380.92*c)) / (2*1380.92)
        if (float(line[0]) == 0 and in_bounds(BV,float(line[4]),const)):
            ngc3680_c.append(BV)
            ngc3680_l.append(float(line[4]))
            lim_ngc3680.append(False)
        elif (float(line[0]) == 1 and in_bounds(BV,float(line[3]),const)):
            ngc3680_c.append(BV)
            ngc3680_l.append(float(line[3]))
            lim_ngc3680.append(True)
    
    ngc3680_c,ngc3680_l = np.array(ngc3680_c),np.log10(np.array(ngc3680_l))
    bv_li.append([ngc3680_c,ngc3680_l])
    upper_lim.append(lim_ngc3680)
    """

    bp_c = []
    bp_l = []
    lim_bp = []
    t = ascii.read('data/beta_pic.txt', delimiter='\t')
    for line in t:
        if in_bounds(float(line[4]),float(line[5]),const):
            bp_c.append(float(line[4]))
            bp_l.append(float(line[5]))
            lim_bp.append(False)
    bp_c,bp_l = np.array(bp_c),np.log10(np.array(bp_l))
    bv_li.append([bp_c,bp_l])
    upper_lim.append(lim_bp)

    fits = get_li_fits(bv_li,upper_lim)

    if (saveToFile):
        pickle.dump(bv_li,open('data/bv_li.p','wb'))
        pickle.dump(upper_lim,open('data/upper_lim.p','wb'))
        pickle.dump(fits,open('data/li_fits.p','wb'))

    return bv_li, upper_lim, fits

        
    
if __name__ == "__main__":
    read_calcium(False,saveToFile=True)
    read_lithium(False,saveToFile=True)

