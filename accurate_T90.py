#### t90.py By Varun/Sujaye
#### Input style given is modified by Vidushi 
#### Issue: print "Total number of photons: ", np.sum(detrend_lc[clipmask]): It takes care of detrend but doesn't subtract mean bkg from grb region
## e.x.	%run accurate_T90.py AS1CZT_GRB171010A_quad_clean.evt GRB171010A --tmark 245358051.58 --tran_start -10.0 --tran_end 110.0
from astropy.io import fits
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clipped_stats, sigma_clip
from matplotlib.backends.backend_pdf import PdfPages
import argparse 
#------------------------------------------------------------------------
# Functions
def trend(time, a, b, c):
    """
    Given time series and some parameters, return the evaluated trend.
    """
    return a * time**2 + b * time + c

#------------------------------------------------------------------------
# Main code

def measure_t90(filename, plotfile, tmin, tmax, usequads, binsize):
    theplot = PdfPages(plotfile)
    quadnames = "".join([ ["A", "B", "C", "D"][quad] for quad in usequads ])
    suptitle = "Quad: {quadnames}, {binsize:0.2f}s binning, File:{filename}".format(quadnames=quadnames, binsize=binsize, filename=filename)

    d = None
    for quad in usequads:
        if d is None:
            d = fits.getdata(filename, quad+1)
        else:
            d = np.hstack((d, fits.getdata(filename, quad+1)))
    #d1 = fits.getdata(filename, 1)
    #d2 = fits.getdata(filename, 1)
    #d = np.hstack((d1, d2))
    #d = np.hstack((d1, d2, d3, d4))
    #d = np.copy(d2)

    d.sort(order='Time')
    select = (d['Time'] >= tmin) & (d['Time'] <= tmax)
    d = d[select]

    #tmin = min([min(d1['time']), min(d2['time'])])
    #tmax = max([max(d1['time']), max(d2['time'])])
    #tmin = np.min(d['Time'])
    #tmax = np.max(d['Time'])

    tbins = np.arange(tmin, tmax+binsize, binsize)
    tmid = 0.5 * (tbins[1:] + tbins[:-1])
    plotx = tmid - tmin
    mask = (tmid < tran_start) | (tmid > tran_end) # True for data points to include
    tranmask = np.copy(mask)

    plt.figure()
    lcurve, tbins = np.histogram(d['Time'], bins=tbins)
    plt.plot(tmid, lcurve, label="Data")

    p0 = [0, 0, np.median(lcurve[mask])]
    x = plotx[mask]
    popt, pcov = curve_fit(trend, x, lcurve[mask], p0)
    first_trend = trend(plotx, *popt)
    plt.plot(tmid, first_trend, label="Initial trend")

    ratio = lcurve/first_trend # Exclude GRB while calculating ratio for sigma-clipping anyway
    rat_clipped = sigma_clip(ratio[mask], sigma=3)
    # rat_clipped.mask is False for data points to include in final fit
    clipmask = np.repeat(True, len(tmid))
    clipmask[mask] = ~rat_clipped.mask
    mask[mask] = ~rat_clipped.mask

    p0 = [0, 0, np.median(lcurve[mask])]
    x = plotx[mask]
    popt, pcov = curve_fit(trend, x, lcurve[mask], p0)
    final_trend = trend(plotx, *popt)
    plt.plot(tmid, final_trend, label='Trend')
    plt.scatter(tmid[~tranmask], lcurve[~tranmask], marker='*', color='k', label='Transient')
    plt.scatter(tmid[~clipmask], lcurve[~clipmask], marker='o', color='r', label='Outliers')
    plt.legend(loc="upper left")
    plt.title("Actual data")
    plt.suptitle(suptitle)
    plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
    plt.xlim((tmin, tmax))
    theplot.savefig()

    plt.figure()
    plt.plot(tmid, ratio)
    plt.title("Ratio")
    plt.suptitle(suptitle)
    plt.scatter(tmid[~tranmask], ratio[~tranmask], marker='*', color='k', label='Transient')
    plt.scatter(tmid[~clipmask], ratio[~clipmask], marker='o', color='r', label='Outliers')
    plt.axhline(1, linestyle='dashed', color='k')
    plt.legend(loc="upper left")
    plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
    plt.xlim((tmin, tmax))
    theplot.savefig()

    plt.figure()
    detrend_lc = lcurve - final_trend
    plt.plot(tmid, detrend_lc)
    plt.title("Detrended LC")
    plt.suptitle(suptitle)
    plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
    plt.axhline(0, linestyle='dashed', color='k')
    plt.xlim((tmin, tmax))
    theplot.savefig()

    plt.figure()
    plt.title(r"T$_{90}$ measurement")
    plt.suptitle(suptitle)
    cs0 = np.cumsum(detrend_lc) / np.sum(detrend_lc)
    cs = np.cumsum(detrend_lc[clipmask]) / np.sum(detrend_lc[clipmask])

    if True:
        # Reject outliers in T90 calculation
        print "Rejecting sigma clip outliers"
        print "Total number of photons: ", np.sum(detrend_lc[clipmask])
        tclip = tmid[clipmask]
        #selclip = (tclip >= tran_start) & (tclip <= tran_end)
        selclip = tclip > 0
        t_trans = tclip[selclip]
        cs_trans = cs[selclip]
    else:
        # No outlier rejection
        print "Keeping sigma clip outliers"
        t_trans = tmid[~tranmask]
        cs_trans = cs[~tranmask]

    # Measure T90 in the transient region
    #cs_check = np.array([0.05, 0.95])
    #t_check = np.interp(cs_check, cs_trans, t_trans) ## Interpolation fails with mutliple intersections
    index_05 = np.max(np.where(cs_trans <= 0.05))
    t05 = np.interp(0.05, cs_trans[index_05:index_05+2], t_trans[index_05:index_05+2]) # Last element of range is not included
    index_95 = np.min(np.where(cs_trans >= 0.95))
    t95 = np.interp(0.95, cs_trans[index_95-1:index_95+1], t_trans[index_95-1:index_95+1]) # Last element of range is not included
    t90 = t95 - t05
    print "T90 is {t90:0.1f} for binsize {binsize:0.2f}".format(t90=t90, binsize=binsize)
    t90text = r"T$_{{90}}$ = {t90:0.2f}".format(t90=t90)

    plt.plot(tmid, cs0, label="Without outlier rejection")
    plt.plot(tmid[clipmask], cs, label="With outlier rejection")
    plt.axhline(0.05, color='gray', alpha=0.5)
    plt.axhline(0.95, color='gray', alpha=0.5)
    plt.axhline(0.00, color='k')
    plt.axhline(1.00, color='k')
    plt.axvline(t05, color='r', label=t90text)
    plt.axvline(t95, color='r')
    #plt.scatter(t_trans, cs_trans)
    plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
    plt.legend(loc="best")
    plt.xlim((tmin, tmax))
    theplot.savefig()

    theplot.close()
    print "Plots saved to ", plotfile

# Input Informations:-----------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("inputfile_single_event",type=str,help='Enter the path of quad clean evt file')
parser.add_argument("GRBNAME", type=str, help='enter the GRB NAME')
parser.add_argument("--tmark",type=float,help='Trigger time for Veto lightcurve')
parser.add_argument("--tran_start",type=float,help='GRB START since trigger time')
parser.add_argument("--tran_end",type=float,help='GRB STOP sincr trigger time')
parser.add_argument("--tbin",type=float,help='Binning time for Veto lightcurve, default=1.0',default=1.0)
args = parser.parse_args()
args.tmark
filename = args.inputfile_single_event
tmin = int(args.tmark) - 200.0*args.tbin
tmax = int(args.tmark) + 300.0*args.tbin

tran_start = args.tmark + args.tran_start
tran_end   = args.tmark + args.tran_end
usequads = [0, 1, 2, 3]
GRBNAME=args.GRBNAME
for binsize in [0.50*args.tbin, 1.00*args.tbin, 1.50*args.tbin, 2.00*args.tbin]:
    plotfile = str(GRBNAME)+"_t90_ABCD_bin{bins:0.3f}.pdf".format(bins=binsize)
    measure_t90(filename, plotfile, tmin, tmax, usequads, binsize)

