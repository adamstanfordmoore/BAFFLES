#!/usr/bin/env python3

"""
Adam Stanford-Moore
5/22/20
Imports the calcium and lithium data and functions.
Lithium data imported int three lists:
bv_li has an array for each cluster of [np.arra([B-V]),np.array([LogEW])]
bv_li = [[[cluster1 B-V],[cluster1 LogEW]],[[cluster2 B-V],[cluster2 LogEW]],....]
upper_lim has an array for each cluster with a boolean for each data point indicating if point is an upper limit
upper_lim [[T,F,T,F,F,F],[F,F,F,F,F]....]
fits stored as [median function,scatter function] for each cluster
"""

from astropy.table import Table
from astropy.io import ascii
from astropy.utils.exceptions import AstropyWarning
import warnings
import numpy as np
import pickle
import fitting as my_fits
import utils
import copy
from os.path import join


def abdor(annotate=False):
    import li_constants as const
    abdor_c = []
    abdor_l = []
    abdor_bverr = []
    abdor_lerr = []
    t = ascii.read(join('data','ab_dor_updated_err_bv_bib.csv'),delimiter=',')
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
        else:
            t[i+1,14] = 'OOR'
            
    if annotate:
        np.savetxt(join("data","ab_dor_updated_err_bv_annotated.csv"),t,delimiter=',',fmt='%s')
    
    
    abdor_c,abdor_l = np.array(abdor_c),np.log10(np.array(abdor_l))
    return abdor_c, abdor_l,abdor_lerr

def tuchor(annotate=False):
    import li_constants as const
    tuchor_c = []
    tuchor_l = []
    tuchor_bverr,tuchor_lerr = [],[]
    t = ascii.read(join('data','tuchor_updated_err_bv_bib.csv'),delimiter=',')
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
        else:
            t[i+1,14] = 'OOR'
    
    if annotate:
        np.savetxt(join("data","tuchor_updated_err_bv_annotated.csv"),t,delimiter=',',fmt='%s')
    
    tuchor_c,tuchor_l = np.array(tuchor_c),np.log10(np.array(tuchor_l))
    return tuchor_c,tuchor_l,tuchor_lerr


#reads in data and generates fits
def read_calcium(fromFile=True,saveToFile=False,fit_degree=0):
    if (fromFile):
        bv_rhk = pickle.load(open(join('data','bv_rhk.p'),'rb'))
        fits = pickle.load(open(join('data','ca_fits.p'),'rb'))
        return bv_rhk,fits
    import ca_constants as const
    warnings.simplefilter('ignore', category=AstropyWarning)
    t = Table.read(join('data','mamajek_table_5.fits'))
    fits = []
    bv_rhk = []
    cluster_index = const.CLUSTER_INDEX
    for i in range(len(cluster_index)):
        c,r = [],[]
        names = []
        if (len(cluster_index[i]) > 2):
            for j in range(0,len(cluster_index[i]),2):
                c += list(t["__B-V_0"][cluster_index[i][j]:cluster_index[i][j+1]])
                r += list(t["logR_HK"][cluster_index[i][j]:cluster_index[i][j+1]])
                names = np.array(t["Name"][cluster_index[i][j]:cluster_index[i][j+1]])
            c = np.array(c).astype(np.float)
            r = np.array(r).astype(np.float)
        else:
            c = np.array(t["__B-V_0"][cluster_index[i][0]:cluster_index[i][1]]).astype(np.float)
            r = np.array(t["logR_HK"][cluster_index[i][0]:cluster_index[i][1]]).astype(np.float)
            names = np.array(t["Name"][cluster_index[i][0]:cluster_index[i][1]])
        
        # omit stars out of bv/metal range
        mask = (const.BV_RANGE[0] <= c) & (c <= const.BV_RANGE[1]) & \
                (const.METAL_RANGE[0] <= r) & (r <= const.METAL_RANGE[1])
        c,r,names = c[mask],r[mask],names[mask]
        
        # average together same stars
        c2,r2 = [],[]
        for name in set(names):
            mask = names == name
            c2.append(np.mean(c[mask]))
            r2.append(np.mean(r[mask]))
        c,r = c2,r2

        if i == const.CLUSTER_NAMES.index(r"$\beta$ Pic"):
            c += [.49,.59,0.803]
            r += [-4.41,-4.37,-4.159]
            #From [Wright 2004, Wright 2004, Gray 2006]
        if i == const.CLUSTER_NAMES.index("Tuc/Hor"):
            c +=  [.52,.6]
            r += [-4.38,-4.33]
            # From [Jenkins 2006,Henry 1996]

        bv_rhk.append([c,r])
        if fit_degree > 0:
            fits.append(my_fits.poly_fit(c,r,n=fit_degree,scatter=True))
        else:
            rhk_fit = np.poly1d([np.median(r)])
            sig_fit = np.poly1d(np.std(my_fits.residuals(c,r,rhk_fit)))
            fits.append([rhk_fit,sig_fit])
    if (saveToFile):
        pickle.dump(bv_rhk,open(join('data','bv_rhk.p'),'wb'))
        pickle.dump(fits,open(join('data','ca_fits.p'),'wb'))
    return bv_rhk,fits

def get_li_fits(bv_li,upper_lim_all):
    import li_constants as const
    fits = []
    for i in range(len(bv_li)):
        fit = None
        if (const.CLUSTER_NAMES[i] == 'Hyades'):
            fit = my_fits.li_dip_fit(bv_li[i][0],bv_li[i][1],upper_lim_all[i],\
                  dip_bv_range=[.398,.513],dip_li_range=[.5,1.6],dip_fit_range=[.39,.52],\
                  edge_box = [.39,.52,1.9,1.95])
        else: 
            poly_order = 2
            fit = my_fits.poly_fit(bv_li[i][0],bv_li[i][1],poly_order,upper_lim_all[i])        
        fits.append(fit)
    fits = my_fits.cluster_scatter_from_stars(bv_li,fits)
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

def make_picklable(fits):
    const = utils.init_constants('lithium')
    fits = copy.deepcopy(fits)
    for c,i in [(c,i) for c in range(len(fits)) for i in range(2)]:
        if type(fits[c][i]) != type(np.poly1d([1])):
            fits[c][i] = fits[c][i](const.BV)
    return fits

def undo_picklable(fits):
    const = utils.init_constants('lithium')
    for c,i in [(c,i) for c in range(len(fits)) for i in range(2)]:
        if type(fits[c][i]) != type(np.poly1d([1])):
            fits[c][i] = my_fits.piecewise(const.BV,fits[c][i])
    return fits

    
def schkolnik_betaPic(annotate=False):
    import li_constants as li_const
    t = np.genfromtxt(join("data","Shkolnik_2017_betaPic_bib.csv"),delimiter=',',dtype=str)
    
    b,l,ul,names = [],[],[],[]
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
    
    if annotate:
        np.savetxt(join("data","Shkolnik_2017_betaPic_annotated.csv"),t,delimiter=',',fmt='%s')
    return b,l,ul,names

def mentuch2008_betaPic(annotate=False):
    import li_constants as li_const
    bp_c = []
    bp_l = []
    lim_bp = []
    names = []
    bp_err = []
    t = np.genfromtxt(join('data','beta_pic_updated_err_2MASS_bib.txt'), delimiter='\t',dtype=str,skip_header=2)
    for i,line in enumerate(t):
        if line[11] != '': continue
        if in_bounds(float(line[4]),float(line[5]),li_const):
            bp_c.append(float(line[4]))
            bp_l.append(float(line[5]))
            lim_bp.append(False)
            names.append(line[10])
            bp_err.append(float(line[6]))
        else:
            t[i,11] = 'OOR'
    bp_c,bp_l = np.array(bp_c),np.log10(np.array(bp_l))
    
    if annotate:
        np.savetxt(join("data","beta_pic_updated_err_2MASS_annotated.txt"),t,delimiter='\t',fmt='%s')
    return bp_c,bp_l,lim_bp,bp_err,names

def merged_betaPic():
    bp_c,bp_l,lim_bp,bp_err,names = mentuch2008_betaPic()
    bp_c2,bp_l2,lim_bp2,names2 = schkolnik_betaPic()
    combined_err = [None]*len(names2)
    for i in range(len(names)):
        if names[i] in names2:
            j = names2.index(names[i])
            #print(names[i],bp_c[i],bp_c2[j],bp_l[i],bp_l2[j],lim_bp[i],lim_bp2[j])
            bv = (bp_c[i]+bp_c2[j])/2
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
    return bp_c2,bp_l2,lim_bp2,combined_err,names2

def alpha_per_lithium():
    import li_constants as const
    c = []
    l = []
    teff_to_bv = my_fits.magic_table_convert('teff','bv')
    t = np.genfromtxt(join("data","alpha_per_balachandra.csv"),delimiter=',',dtype=str)
    for line in t:
        if line[11] != '': continue
        bv = teff_to_bv(float(line[1]))
        ew = float(line[7])
        if in_bounds(bv,ew,const):
            c.append(bv)
            l.append(ew)
    c,l = np.array(c),np.log10(l)
    return c,l,[False]*len(c)

def read_lithium(fromFile=True,saveToFile=False):
    if (fromFile):
        bv_li = pickle.load(open(join('data','bv_li_all.p'),'rb'))
        upper_lim = pickle.load(open(join('data','upper_lim_all.p'),'rb'))
        fits = pickle.load(open(join('data','li_fits_all.p'),'rb'))
        fits = undo_picklable(fits)
        return bv_li,upper_lim,fits
    
    import li_constants as const
    warnings.simplefilter('ignore', category=AstropyWarning)
    bv_li = []
    upper_lim = []
 

    t = ascii.read(join('data','ngc2264_lithium_bv.csv'),delimiter=';')
    ngc2264_c = []
    ngc2264_l = []
    for line in t[2:]:  
        if (float(line[6]) != -999.0 and in_bounds(float(line[7]),1000*float(line[6]),const)):
            ngc2264_c.append(float(line[7]))
            ngc2264_l.append(1000*float(line[6]))

    t = ascii.read(join('data','ngc2264_lithium2.txt'), delimiter=',')
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
    
    bp_c,bp_l,lim_bp,_,_ = merged_betaPic()
    bv_li.append([bp_c,bp_l])
    upper_lim.append(lim_bp)

   
    ic2602_c = []
    ic2602_l = []
    t = ascii.read(join('data','ic2602_lithium.txt'), delimiter=',')
    for line in t:
            if in_bounds(line[1],line[3],const):
                    ic2602_c.append(line[1])
                    ic2602_l.append(line[3])
    ic2602_c,ic2602_l = np.array(ic2602_c),np.log10(np.array(ic2602_l))
    bv_li.append([ic2602_c,ic2602_l])
    upper_lim.append([False]*len(ic2602_c))
    

    aper_c,aper_l,lim_aper = alpha_per_lithium()
    bv_li.append([aper_c,aper_l])
    upper_lim.append(lim_aper)
   
    pleiades_c = []
    pleiades_l = []
    lim_p = []
    t = ascii.read(join('data','pleiades_lithium.tsv'), delimiter=';')
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
    
    m35_c = []
    m35_l = []
    lim_m35 = []
    t = ascii.read(join('data','M35_data.txt'))
    for line in t:
        if in_bounds(float(line[1]),float(line[2]),const):
            m35_c.append(float(line[1]))
            m35_l.append(float(line[2]))
            lim_m35.append(line[-1]!=0)
    m35_c,m35_l = np.array(m35_c),np.log10(np.array(m35_l))
    bv_li.append([m35_c,m35_l])
    upper_lim.append(lim_m35)
    
    m34_c = []
    m34_l = []
    lim_m34 = []
    t = ascii.read(join('data','m34_lithium.txt'), delimiter=',')
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

    coma_c = []
    coma_l = []
    lim_coma = []
    t = ascii.read(join('data','coma_berenices.txt'), delimiter=',')
    for line in t:
        if line[-1] != '': continue
        if in_bounds(float(line[2]),float(line[4]),const):
            coma_c.append(float(line[2]))
            coma_l.append(float(line[4]))
            lim_coma.append(line[7]==1)
    coma_c,coma_l = np.array(coma_c),np.log10(np.array(coma_l))
    bv_li.append([coma_c,coma_l])
    upper_lim.append(lim_coma)


    hyades_c = []
    hyades_l = []
    lim_h = []
    t = ascii.read(join('data','hyades_lithium.tsv'), delimiter=';')
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
    t = ascii.read(join('data','m67_lithium_eric_edits.txt'), delimiter=',')
    for line in t:
        # filter from Eric's notes
        if line[-1] != 'single member' and line[-1] != 'binary member':
            continue
        if line[-1] == 'binary member' and float(line[4]) < 0:
            continue
        
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

    """
    uma_c = []
    uma_l = []
    lim_uma = []
    t = ascii.read(join('data','UMa.csv'), delimiter=',')
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
    t = ascii.read(join('data','ngc_3680.csv'),delimiter=',')
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

    fits = get_li_fits(bv_li,upper_lim)

    if (saveToFile):
        pickle.dump(bv_li,open(join('data','bv_li_all.p'),'wb'))
        pickle.dump(upper_lim,open(join('data','upper_lim_all.p'),'wb'))
        pickle.dump(make_picklable(fits),open(join('data','li_fits_all.p'),'wb'))

    return bv_li, upper_lim, fits

        
    
if __name__ == "__main__":
    read_calcium(False,saveToFile=True)
    read_lithium(False,saveToFile=True)

