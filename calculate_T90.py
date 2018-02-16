##=================================================================================================================================
## Calculate T90
## By Vidushi Sharma, 17 Dec 2017
##=================================================================================================================================
## Step 1: Banana pixel sorted out and livetime corrected
## e.x.: %run calculate_T90.py AS1CZT_GRB171010A_quad_clean.evt AS1cztbadpix20160908v01.fits AS1C03_016T01_9000001596_11010cztM0_level2_quad_livetime.fits GRBname --tmark 245358051.58 
##=================================================================================================================================

import numpy as np
import os,subprocess,argparse 
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.io import fits
import math
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages

# Input Informations:-----------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("inputfile_single_event",type=str,help='Enter the path of quad clean evt file')
parser.add_argument("inputfile_banana_pix",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("inputfile_livetime",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("outfile", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("GRBNAME", type=str, help='enter the GRB NAME')
parser.add_argument("--tmark",type=float,help='Trigger time for lightcurve')
parser.add_argument("--tbin",type=float,help='Binning time for lightcurve, default=1.0',default=1.0)
parser.add_argument("-p", "--plottype", help="Type of output file to produce (png, pdf)", type=str, default='pdf')

args = parser.parse_args()
args.tmark
print args.tbin
tmin = int(args.tmark) - 300.0*args.tbin
tmax = int(args.tmark) + 500.0*args.tbin
print(tmin,tmax,args.tmark)
print ('\n \n Required Input files: \n (1)	*_quad_clean.dblevt \n (2)	Caldb badpix file, AS1cztbadpix20160908v01.fits \n (3)	Livetime fits file according tbin or type of GRB	\n')
print('trigger time = %.2f, tmin = %.2f, tmax = %.2f	\n' %(args.tmark, tmin, tmax))
tbins = np.arange(tmin, tmax + args.tbin, args.tbin)

evtfile	=	fits.open(args.inputfile_single_event)
badfile	=	fits.open(args.inputfile_banana_pix)
livefile	=	fits.open(args.inputfile_livetime)

hist	=	[0.0]*int((tmax-tmin)/args.tbin)
for Q_n in range(1,5):
# Open and read the Input file: single event file
	print("		For Quadrant %d:" % Q_n)
	evtable	=	evtfile[Q_n].data
	evtime	=	evtable['Time']
	evtdet	=	evtable['DetID']
	evtpix	=	evtable['pixID']
	evtene	=	evtable['ENERGY']
	Vtable	=	Table([evtime,evtdet,evtpix,evtene], names=('Time','DetID','pixID','Energy'))
	start_match 	=	np.where(evtime >= tmin)[0]
	start	=	start_match[0]
	stop_match	=	np.where(evtime >= tmax)[0]
	stop	= stop_match[0]

	evtable	=	Vtable[start:stop]
	evtime	=	evtable['Time']
	evtdet	=	evtable['DetID']
	evtpix	=	evtable['pixID']
	energy	=	evtable['Energy']
	print("Event length = %d" % len(evtime))
	"""
	badtable=	badfile[Q_n].data	
	badind	=	np.where(badtable['PIX_FLAG']==1)
	baddet	=	(badtable['DETID'])[badind[0]]
	badpix	=	(badtable['PIXID'])[badind[0]]
	badflag	=	(badtable['PIX_FLAG'])[badind[0]]
	print('Length of badpix with flag 1 = %d'%len(badind[0]))

	goodevt	=	[]
	badevt	=	[]
	for i in range (0,len(evtdet)):
		flag = 0
		for j in range (0,len(baddet)):
			if((evtdet[i]==baddet[j]) and (evtpix[i]==badpix[j])):
				badevt.append([i,j])
				flag=1
				break
		if(flag==0):
			goodevt.append(i)
			
	badevt	=	np.array(badevt)
	goodevt	=	np.array(goodevt)
	print("banana pixel length = %d" % len(badevt))
	print("cleaned event length = %d" % len(goodevt))
	Event_table	=	evtable[goodevt]"""
	Event_table	=	evtable
	Q1_hist, bin_edges = np.histogram(Event_table['Time'] , bins=tbins)
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0 
	Q1_hist	=	Q1_hist/float(args.tbin)
	#plt.plot(center_time,Q1_hist/float(args.tbin))
	#print(Q1_hist)
## LIVETIME correction is done: For 1 s lightcurve dividing the 1 s collected counts by fracexp of livetime fits.
##===============================================================================================================================
	livetable	=	livefile[Q_n].data
	time	=	livetable['TIME']
	fracexp	=	livetable['FRACEXP']
	ltable	=	Table([time, fracexp], names=('Time', 'Fracexp'))

	lstart_match 	=	np.where(time >= tmin)[0]
	lstart	=	lstart_match[0]
	lstop_match 	=	np.where(time >= tmax)[0]
	lstop	=	lstop_match[0]

	evtlive	=	ltable[lstart:lstop]
	ltime	=	evtlive['Time']
	lfrac	=	evtlive['Fracexp']
	Q1_hist	=	Q1_hist
	#/lfrac
	#print(lfrac)
	#print ('Livetime correction done \n')
	hist = hist + Q1_hist
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0 
############################################################################################
Counts	=	hist
Time	=	center_time
print(Counts[0],Time[0])			
## SELECT PRE BACKGROUND INTERVAL
print("\n SELECT PRE GRB BACKGROUND START TO STOP INTERVAL *********************")
# Simple mouse click function to store coordinates
def onclick(event):
	global ix, iy
	ix, iy = event.xdata, event.ydata
	# assign global variable to access outside of function
	global coords
	coords.append((ix, iy))

 	# Disconnect after 2 clicks
        if len(coords) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(Time, Counts)
ax.set_xlim(tmin, tmax )
coords = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
pre_bkg_min=coords[0][0]
pre_bkg_stop=coords[1][0]	
print("Pre-bkg start = %.2f ,Pre-bkg stop= %.2f" %(pre_bkg_min, pre_bkg_stop))

## SELECT GRB INTERVAL
print("\n SELECT THE GRB START TO STOP INTERVAL ********************************")
def onclick(event):
	global ix1, iy1
	ix1, iy1 = event.xdata, event.ydata
	global coords_grb
	coords_grb.append((ix1, iy1))
	if len(coords_grb) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(Time, Counts)
ax.set_xlim(pre_bkg_min,tmax)
coords_grb = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
grb_start=coords_grb[0][0]
grb_stop=coords_grb[1][0]	
print("GRB start = %.2f , GRB stop = %.2f" %(grb_start,grb_stop))

## SELECT GRB INTERVAL
print("\n SELECT THE POST BKG START TO STOP INTERVAL ********************************")
def onclick(event):
	global ix2, iy2
	ix2, iy2 = event.xdata, event.ydata
	global coords_2
	coords_2.append((ix2, iy2))
	if len(coords_2) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(Time, Counts)
ax.set_xlim(pre_bkg_min,tmax)
coords_2 = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
post_bkg_start=coords_2[0][0]
post_bkg_stop=coords_2[1][0]	
print("Post-bkg start = %.2f , Post-bkg stop = %.2f" %(post_bkg_start, post_bkg_stop))
print(post_bkg_start, post_bkg_stop)
##############################################################################################

pre_bkg_stamp	=	np.where(center_time >= pre_bkg_min)[0]
t1	=	pre_bkg_stamp[0]
pre_bkg_stamp2	=	np.where(center_time >= pre_bkg_stop)[0]
t2	=	pre_bkg_stamp2[0]
grb_start_stamp	=	np.where(center_time >= grb_start)[0]
tgrb1	=	grb_start_stamp[0]
grb_start_stamp2	=	np.where(center_time >= grb_stop)[0]
tgrb2	=	grb_start_stamp2[0]
post_bkg_stamp	=	np.where(center_time >= post_bkg_start)[0]
t3	=	post_bkg_stamp[0]
post_bkg_stamp2	=	np.where(center_time >= post_bkg_stop)[0]
t4	=	post_bkg_stamp2[0]

############ VEDANT CODE bit modified below #####################################
b_counts=np.hstack((Counts[t1:t2],Counts[t3:t4])) ###  Pre Post Bkg
b_time=np.hstack((Time[t1:t2],Time[t3:t4])) ###  Pre Post Bkg

g_counts=np.hstack((Counts[tgrb1:tgrb2])) ###  GRB
g_time=np.hstack((Time[tgrb1:tgrb2])) ### GRB

pfit = np.polyfit(b_time, b_counts,2)  ## fitting with quadratic
bkg=np.polyval(pfit,Time)
subs_data=Counts-bkg
result=gaussian_filter1d(subs_data,sigma=0, axis=-1, order=0)

peak_count=max(g_counts)-np.average(bkg);
peak_time=g_time[np.where(g_counts==max(g_counts))[0][0]]

t_90_data=result[tgrb1:tgrb2]  ###  Manually add the range for GRB Duration 
t_90_accum=np.add.accumulate(t_90_data)
t_90_time=Time[tgrb1:tgrb2] ###  Manually add the range for GRB Duration
max_acc=max(t_90_accum)

fracout = np.array([max_acc*0.05, max_acc*0.95])
tout = np.interp(fracout, t_90_accum, t_90_time)
t90 = tout[1] - tout[0]

def gen_plot(tmark,peak_count,peak_time,max_acc,t90,Time,Counts,bkg,result,subs_data):
	print ("peak count from background: %f at %f s astrosat seconds" % (peak_count, peak_time))
	print ("mean background count: %f" % np.mean(bkg))
	print ("total_count accumulated: %f " % max_acc)
	print ("T90 Calculated for %s: %f" % (args.GRBNAME,t90))
	if args.plottype == 'pdf':
		plotfile = PdfPages("{stem}_T90.pdf".format(stem=args.outfile))
 
	plt.figure()
	plt.subplot(211)
	plt.title('Count Rate Plot')
	plt.ylabel('Count (Counts/sec)', color='red')
	plt.xlabel('Time (sec)', color='red')
	plt.plot(Time,Counts)
	plt.plot(Time,bkg)
	plt.plot(Time,result)
	plt.plot(Time,subs_data)
	plt.text(tmin,peak_count+0.07*peak_count, "Peak Count: %f at %f Astrosat Second" % (peak_count,peak_time))
	plt.text(tmin,peak_count+0.14*peak_count, "Total Count : %f " % max_acc)
	plt.text(tmin,peak_count+0.21*peak_count, "Mean background Count: %f" % np.mean(bkg)) 
	plt.grid()
	if args.tmark > 0: plt.axvline([args.tmark], color='gray', linestyle='dashed', label='Marked times')
	plt.subplot(212)
	plt.title('Accumulated Count plot to Calculate T90')
	plt.tight_layout()
	plt.plot(t_90_time,t_90_accum)
	plt.ylabel('Accumulated Count (Counts/sec)', color='red')
	plt.xlabel('Time (sec)', color='red')
	plt.text(t_90_time[0],max_acc-0.1*max_acc,"T90 Calculated for %s: %f" % (args.GRBNAME,t90))   
	#plt.show()
	#GRBNAME=args.GRBNAME
	#GRBNAME=str(GRBNAME)+str('_T90.pdf')  ######## It doesnt save anything
	#print (GRBNAME)
	#return plt.savefig(GRBNAME), plt.show()
	#plt.show()
	if args.plottype == 'pdf':
		plotfile.savefig()
	else:
		plt.show()	
		plt.savefig(args.outfile + "_T90." + args.plottype)
	#GRBNAME=args.GRBNAME
	#print(GRBNAME)
	#GRBNAME=str('AS1CZT_')+GRBNAME+str('_VetoLightCurve.png')
	#print(GRBNAME)
	#GRBLOCATION=args.GRBLOCATION
	#print(GRBLOCATION)
	#GRBNAME=str(GRBLOCATION)+str(GRBNAME)
	#print(GRBNAME)
	#plt.savefig(GRBNAME)
	plt.show()	
	if args.plottype == 'pdf':
    		plotfile.close()



gen_plot(args.tmark,peak_count,peak_time,max_acc,t90,Time,Counts,bkg,result,subs_data)

## rsync -r GRB160422A astrosat@purva:/wusr3/project/astrosat/public_html/czti_grb




