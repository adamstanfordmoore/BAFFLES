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
import utils
import copy

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
        
        # omit stars out of bv/metal range
        mask = (const.BV_RANGE[0] <= c) & (c <= const.BV_RANGE[1]) & \
                (const.METAL_RANGE[0] <= r) & (r <= const.METAL_RANGE[1])
        c,r = c[mask],r[mask]
        
        bv_rhk.append([c,r])
        #fits.append(my_fits.constant_fit(r))
        #fits.append([np.poly1d([np.median(r)]),np.poly1d([0.1])])
        #fits.append(my_fits.poly_fit(c,r,n=fit_degree,scatter=True))
        rhk_fit = np.poly1d([np.median(r)])
        sig_fit = np.poly1d(np.std(my_fits.residuals(c,r,rhk_fit)))
        fits.append([rhk_fit,sig_fit])
    if (saveToFile):
        pickle.dump(bv_rhk,open('data/bv_rhk.p','wb'))
        pickle.dump(fits,open('data/ca_fits.p','wb'))
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
        
        #if (const.CLUSTER_NAMES[i] in ['UMa','NGC3680']): #these cluster need linear fits
        #    poly_order = 1
        #fit = my_fits.pwise_fit(bv_li[i][0],bv_li[i][1],upper_lim=upper_lim_all[i],guess_fit=my_fits.poly_fit(bv_li[i][0],bv_li[i][1]),x_method='bin')
        
        
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
    #i = const.CLUSTER_NAMES.index('Hyades')
    #fits = copy.deepcopy(fits)
    for c,i in [(c,i) for c in range(len(fits)) for i in range(2)]:
        if type(fits[c][i]) != type(np.poly1d([1])):
            fits[c][i] = fits[c][i](const.BV)
    return fits

def undo_picklable(fits):
    const = utils.init_constants('lithium')
    #i = const.CLUSTER_NAMES.index('Hyades')
    for c,i in [(c,i) for c in range(len(fits)) for i in range(2)]:
        if type(fits[c][i]) != type(np.poly1d([1])):
            fits[c][i] = my_fits.piecewise(const.BV,fits[c][i])
    #fits[i][0] = my_fits.piecewise(const.BV,fits[i][0])
    #fits[i][1] = my_fits.piecewise(const.BV,fits[i][1])

def read_lithium(fromFile=True,saveToFile=False):
    if (fromFile):
        bv_li = pickle.load(open('data/bv_li_all.p','rb'))
        upper_lim = pickle.load(open('data/upper_lim_all.p','rb'))
        fits = pickle.load(open('data/li_fits_all.p','rb'))
        undo_picklable(fits)
        #bv_li = pickle.load(open('data/bv_li_gaia.p','rb'))
        #upper_lim = pickle.load(open('data/upper_lim_gaia.p','rb'))
        #fits = pickle.load(open('data/li_fits_gaia.p','rb'))
        return bv_li,upper_lim,fits
    
    import li_constants as const
    bv_li = []
    upper_lim = []
 
    """
    teff_to_bv = my_fits.magic_table_convert('teff','bv')
    #reads in easy to read cluster table
    def readCluster(filePath,bv_index,EW_index,teff=False,UL_index=None,delimeter=';',mem_indices=[-1,-1],prob_threshold=.9):
        c,l,ul = [],[],[]
        t = ascii.read(filePath, fill_values='blank',delimiter=delimeter)
        rejected,kept = 0,0
        for line in t[2:]:
            #check membership
            if mem_indices[0] != -1 and (not utils.isFloat(line[mem_indices[0]]) \
                or float(line[mem_indices[0]]) < prob_threshold): 
                #print("Rejecting",line[mem_indices[0]])
                rejected += 1
                continue
            if mem_indices[1] != -1 and (not utils.isFloat(line[mem_indices[1]]) \
                or float(line[mem_indices[1]]) < prob_threshold): 
                #print("B-Rejecting",line[mem_indices[1]])
                rejected += 1
                continue
            if not utils.isFloat(line[bv_index]) or not utils.isFloat(line[EW_index]):
                #print("not float",line[bv_index],line[EW_index])
                rejected += 1
                continue
            bv = teff_to_bv(float(line[bv_index])) if teff else float(line[bv_index]) 
            ew = float(line[EW_index])
            if not in_bounds(bv,ew,const): 
                #print "not in bounds",bv,np.log10(ew)
                rejected += 1
                continue
            kept += 1
            c.append(bv)
            l.append(ew) 
            if UL_index and line[UL_index] == '<': #upper_lim
                ul.append(True)
            else:
                ul.append(False)
        print "rejected",rejected,"Kept",kept
        c,l = np.array(c),np.log10(np.array(l))
        bv_li.append([c,l])
        upper_lim.append(ul)

    readCluster('data/cha_I_lithium_gaia.txt',bv_index=3,EW_index=5,teff=True,UL_index=None)
    readCluster('data/ic4665_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[17,-1])
    readCluster('data/ic2602_lithium_gaia.txt',bv_index=5,EW_index=12,teff=True,UL_index=11,mem_indices=[19,-1])
    #readCluster('data/ngc2451_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[17,-1])
    readCluster('data/ngc2451_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[-1,18])
    readCluster('data/ngc2547_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[17,-1])
    readCluster('data/ic2391_lithium_gaia.txt',bv_index=5,EW_index=12,teff=True,UL_index=11,mem_indices=[19,-1])
    readCluster('data/ngc2516_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[17,-1])
    readCluster('data/ngc6633_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9,mem_indices=[17,-1],prob_threshold=0.8)



    fits = get_li_fits(bv_li,upper_lim)

    if (saveToFile):
        pickle.dump(bv_li,open('data/bv_li_gaia.p','wb'))
        pickle.dump(upper_lim,open('data/upper_lim_gaia.p','wb'))
        pickle.dump(fits,open('data/li_fits_gaia.p','wb'))

    return bv_li, upper_lim, fits
    """
  





    """
    t = ascii.read('data/NGC2264_bouvier.txt',delimiter=';')
    ngc2264_c = []
    ngc2264_l = []
    print t
    print t[0]
    """
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
    t = ascii.read('data/m67_lithium_eric_edits.txt', delimiter=',')
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
    """
    teff_to_bv = my_fits.magic_table_convert('teff','bv')
    #reads in easy to read cluster table
    def readCluster(filePath,bv_index,EW_index,teff=False,UL_index=None,delimeter=';'):
        c,l,ul = [],[],[]
        t = ascii.read(filePath, delimiter=delimeter)
        for line in t[2:]:
            bv = teff_to_bv(float(line[bv_index])) if teff else float(line[bv_index]) 
            ew = float(line[EW_index])
            if not in_bounds(bv,ew,const): continue
            c.append(bv)
            l.append(ew) 
            if UL_index and line[UL_index] == '<': #upper_lim
                ul.append(True)
            else:
                ul.append(False)
        c,l = np.array(c),np.log10(np.array(l))
        bv_li.append([c,l])
        upper_lim.append(ul)

    readCluster('data/cha_I_lithium_gaia.txt',bv_index=3,EW_index=5,teff=True,UL_index=None)
    readCluster('data/ic4665_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9)
    readCluster('data/ic2602_lithium_gaia.txt',bv_index=5,EW_index=12,teff=True,UL_index=11)
    readCluster('data/ngc2451_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9)
    readCluster('data/ngc2547_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9)
    readCluster('data/ic2391_lithium_gaia.txt',bv_index=5,EW_index=12,teff=True,UL_index=11)
    readCluster('data/ngc2516_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9)
    readCluster('data/ngc6633_lithium_gaia.txt',bv_index=3,EW_index=10,teff=True,UL_index=9)
    """


    fits = get_li_fits(bv_li,upper_lim)

    if (saveToFile):
        pickle.dump(bv_li,open('data/bv_li_all.p','wb'))
        pickle.dump(upper_lim,open('data/upper_lim_all.p','wb'))
        pickle.dump(make_picklable(fits),open('data/li_fits_all.p','wb'))
        #pickle.dump(bv_li,open('data/bv_li_gaia.p','wb'))
        #pickle.dump(upper_lim,open('data/upper_lim_gaia.p','wb'))
        #pickle.dump(fits,open('data/li_fits_gaia.p','wb'))

    return bv_li, upper_lim, fits

        
    
if __name__ == "__main__":
    read_calcium(False,saveToFile=True)
    read_lithium(False,saveToFile=True)

