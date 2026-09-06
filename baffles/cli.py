"""Command-line interface for BAFFLES.  Installed as the ``baffles`` console script."""
import sys
from baffles import utils
from baffles.core import baffles_age


def main(argv=None):
    if argv is None:
        argv = sys.argv
    const = utils.init_constants('lithium')
    err = "Usage:  baffles -bmv <B-V> -rhk <Log10(R\'HK)> -li <EWLi> [-bmv_err <> -li_err <> -ul -maxAge <13000> -plot -s -savePlot -filename <> -help]"
    
    help_msg = "\n\
    -bmv corrected (B-V)o of the star (optional for calcium)\n\
    -rhk <> log base 10 of the R'HK value \n\
    -li <> EW measure in milli-angstroms (=0.1 pm) \n\
    \noptional flags:\n\
    -ul indicates that the log(EW/mA) reading is an upper-limit reading. \n\
    -maxAge allows user to input max posible age of star (Myr) if upper-limit flag is used. default is %d \n\
    -bmv_err <float uncertainity> provides the uncertainty on B-V with default %.2f \n\
    -li_err <float uncertainity> provides the uncertainty in LiEW measurement with default %dmA \n\
    -s saves posteriors as .csv as age, probability in two 1000 element columns.\n\
    -plot plots and shows the PDF. \n\
    -savePlot saves the plotted posteriors to a pdf. \n\
    -filename <name of file w/o extension> name of files to be saved: name.pdf is graphs, name_calcium.csv/name_lithium.csv/name_product.csv are posterior csv files for calcium/lithium/product respectively.  \n\
    -help prints this message \n" % (const.GALAXY_AGE,const.BV_UNCERTAINTY, const.MEASURE_ERR)
    
    if (len(argv) < 3 or '-help' in argv):
        print(err,help_msg)
        sys.exit()
    bv = None
    bv_err,li_err = None,None
    rhk, li = None,None
    save = False
    showPlots = False
    fileName = 'baffles'
    savePlots = False
    upperLim = False
    maxAge = const.GALAXY_AGE
    valid_flags = ['-bmv','-rhk','-li','-li_err','-bmv_err','-plot','-savePlot','-ul','-maxAge','-s','-filename','-help']
    extra_flags = ['-Plot','-plots','-Plots','-savePlots','-saveplots','-saveplot','-UL','-save']
    for i,ar in enumerate(argv[1:]):
        if ar not in valid_flags and ar not in extra_flags \
            and not utils.isFloat(ar) and argv[i] != '-filename':
            print("Invalid flag '" + ar + "'. Did you mean one of these:")
            print(valid_flags)
            sys.exit()
    try:
        if ('-bmv' in argv):
            bv = float(argv[argv.index('-bmv') + 1])
        if ('-rhk' in argv):
            rhk = float(argv[argv.index('-rhk') + 1])
            from baffles import ca_constants as const
            if bv is not None and (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print("B-V out of range. Must be in range " + str(const.BV_RANGE))
                sys.exit()
            if (not (const.METAL_RANGE[0] <= rhk <= const.METAL_RANGE[1])):
                print("Log(R\'HK) out of range. Must be in range " + str(const.METAL_RANGE))
                sys.exit()
        if ('-li' in argv):
            li = float(argv[argv.index('-li') + 1])
            from baffles import li_constants as const
            if 0 < li < 3:
                print("Interpretting LiEW as log(LiEW)")
                li = 10**li
            
            if (not (const.BV_RANGE[0] <= bv <= const.BV_RANGE[1])):
                print("B-V out of range. Must be in range " + str(const.BV_RANGE))
                sys.exit()
            if (not (const.METAL_RANGE_LIN[0] <= li <= const.METAL_RANGE_LIN[1])):
                print("Li EW out of range. Must be in range " + str(const.METAL_RANGE) + " mA")
                sys.exit()
        if ('-s' in argv or '-save' in argv):
            save = True
        if ('-plot' in argv or '-Plot' in argv or '-plots' in argv or '-Plots' in argv):
            showPlots = True
        if ('-saveplot' in argv or '-saveplots' in argv or '-savePlots' in argv or '-savePlot' in argv):
            savePlots = True
        if ('-filename' in argv):
            fileName = argv[argv.index('-filename') + 1]
        if ('-bmv_err' in argv):
            bv_err = float(argv[argv.index('-bmv_err') + 1])
        if ('-li_err' in argv):
            li_err = float(argv[argv.index('-li_err') + 1])
        if ('-ul' in argv or '-UL' in argv):
            upperLim = True
        if ('-maxAge' in argv):
            maxAge = float(argv[argv.index('-maxAge') + 1])
    except IndexError:
        print(err)
    except ValueError:
        print(err)
    
    baffles_age(bv,rhk,li,bv_err,li_err,upperLim,maxAge,fileName,showPlots=showPlots,
                savePlots=savePlots, savePostAsText=save)


if __name__ == "__main__":
    main()
