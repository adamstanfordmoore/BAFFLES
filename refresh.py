"""
Adam Stanford-Moore
2/11/19
"""
import fitting as my_fits
import baffles
import readData
import utils
import datetime

def main():
    date = datetime.datetime.now().strftime("%m%d%y")
    bv_m,fits = readData.read_calcium(fromFile=False,saveToFile=True)
    _,res_arr = my_fits.get_fit_residuals(bv_m,fits,'calcium',None,li_range=None,linSpace=False,scale_by_std= True)
    my_fits.fit_histogram('calcium',residual_arr=res_arr,fromFile=False,saveToFile=True)
   
    baf = baffles.age_estimator('calcium',default_grids=False)
    baf.make_grids(bv_m,fits,medianSavefile='grids/median_rhk_'+date,\
            sigmaSavefile='grids/sigma_rhk_'+date,setAsDefaults=True)

    const = utils.init_constants('lithium')    
    bv_m, upper_lim, fits = readData.read_lithium(fromFile=False,saveToFile=True)
    my_fits.MIST_primordial_li(ngc2264_fit=fits[const.CLUSTER_NAMES.index('NGC2264')][0],fromFile=False, saveToFile=True)
    _,res_arr= my_fits.get_fit_residuals(bv_m,fits,'lithium',upper_lim,li_range=None,linSpace=False)
    my_fits.fit_histogram('lithium',residual_arr=res_arr,fromFile=False,saveToFile=True)
    
    baf2 = baffles.age_estimator('lithium',default_grids=False)
    baf2.make_grids(bv_m,fits,upper_lim,'grids/median_li_'+date,'grids/sigma_li_'+date,setAsDefaults=True)

if  __name__ == "__main__":
    main()
