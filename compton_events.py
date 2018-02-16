##====================================================================================================================================
## Calculate and plot comp and non comp events 
## Vidushi Sharma, 15th Dec 2017
##Last modified-Priyanka Shahane,10 feb 2018
##=====================================================================================================================================
## Banana pixel sorted out and livetime correction is applied in lightcurve
## e.x.: %run comp_noncomp.py AS1C03_016T01_9000001596_11010cztM0_level2_quad_clean.dblevt AS1cztbadpix20160908v01.fits AS1C03_016T01_9000001596_11010cztM0_level2_quad_livetime.fits --tmark 245358051.58
##====================================================================================================================================

import numpy as np
import argparse
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import math
from matplotlib.backends.backend_pdf import PdfPages

np.seterr(over='ignore')

## Input Informations:------------------------------------------------------------------------------------------------------------------
print ('\n \n Required Input files: \n (1)	*_quad_clean.dblevt \n (2)	Caldb badpix file, AS1cztbadpix20160908v01.fits \n (3)	Livetime fits file according tbin or type of GRB	\n')
## -------------------------------------------------------------------------------------------------------------------------------------
parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_double_event",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("inputfile_banana_pix",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("inputfile_livetime",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("outfile", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("--tmark",type=float,help='Trigger time for Veto lightcurve')
parser.add_argument("--tbin",type=float,help='Binning time for Veto lightcurve, default=1.0',default=1.0)
parser.add_argument("--emin",type=float,help='Lower energy limit, default=125',default=100.0)
parser.add_argument("--emax",type=float,help='Upper energy limit, default=250',default=400.0)
parser.add_argument("-p", "--plottype", help="Type of output file to produce (png, pdf)", type=str, default='pdf')
args	=	parser.parse_args()
args.tmark

tmin = int(args.tmark) - 200.0*args.tbin
tmax = int(args.tmark) + 500.0*args.tbin
print('trigger time = %.2f, tmin = %.2f, tmax = %.2f	\n' %(args.tmark, tmin, tmax))
tbins = np.arange(tmin, tmax + args.tbin, args.tbin)

dblevtfile	=	fits.open(args.inputfile_double_event)
badfile		=	fits.open(args.inputfile_banana_pix)
livefile	=	fits.open(args.inputfile_livetime)
	
Q_hist_sum	=	[0.0]*int((tmax-tmin)/args.tbin)
Q_max	=	[0.0]*1
Time_sum	=	[]
if args.plottype == 'pdf':
    plotfile = PdfPages("{stem}_compton_events.pdf".format(stem=args.outfile))

for Q_n in range(1, 5):
# Open and read the Input file:----------------------------------------------------------------------------------
	print("		For Quadrant %d:" % Q_n)
	dblevt	=	dblevtfile[Q_n].data
	evtime	=	dblevt['Time']
	evtdet	=	dblevt['DetID']
	evtpix	=	dblevt['pixID']
	evtene	=	dblevt['ENERGY']
	evtick	=	dblevt['CZTNTICK']
	evtdetX	=	dblevt['DETX']
	evtdetY	=	dblevt['DETY']
	Vtable	=	Table([evtime,evtdet,evtpix,evtene,evtick,evtdetX,evtdetY],	names=('Time','DetID','pixID','Energy','cztnTick','detX','detY'))

	for i in range(0,len(evtime)):
		if (evtime[i] >= tmin):
			start = i
			break

	for j in range(0,len(evtime)):
		if (evtime[j] >= tmax):
			stop = j
			break

	evtable	=	Vtable[start:stop]
	evtime	=	evtable['Time']
	evtdet	=	evtable['DetID']
	evtpix	=	evtable['pixID']
	print("Event length = %d" % len(evtime))
	"""
	## banana pix correction: Whichever detID and pixel ID has flag 1 in caldb badpix fits file, remove the same detID and pixel ID from the input file
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
	cleaned_evt = evtable[goodevt]"""
	cleaned_evt = evtable
	Time	=	cleaned_evt['Time']
	DetID	=	cleaned_evt['DetID']
	PixID	=	cleaned_evt['pixID']
	Energy	=	cleaned_evt['Energy']
	Tick	=	cleaned_evt['cztnTick']
	DetX	=	cleaned_evt['detX']
	DetY	=	cleaned_evt['detY']
	
####=====================================================================================================================================
####	COMPTON CRITERIA	: 
#### (1)	Finding the adjacent pixels with condition that distance between 2 pixel, 1 cm	<=	distance	<=	1.8 cm
####	Reason: pixel size 2.5 cm, if pixel center point is (0,0) then max next pixel edge is (1.25,1.25) diagonally and min next pixel edge ####		is (0,1.25) or (1.25,0). The max distance is 1.76 cm and min distance is 1.25.
#### (2)	1	<	absorbed energy(h) / scattered energy(l)	<	6 
####	Reason: Scattering pixel will have low energy, whereas absorbing pixel will absorb higher energy photon. 
####	If we want to consider photons with min energy 50 keV and max energy 300 keV. Also if compton effect has high cross-section ##	####	for upcoming photon(total energy) in 100 keV to 400 keV. If low energy(scattered) photon at 50 keV, we detect photon in next pixel 		with min energy and max energy 50 and 300 keV respectively, then 1 <= h/l <= 6
######## WHY NOT EMAX = 500 keV ????
	dist	=	[]
	match_1st	=	[]
	match_2nd	=	[]	
	for i in range(0,	len(cleaned_evt)-1):
	## Selecting 2 same tick events and note distance between them
		if(	(Tick[i] - Tick[i+1]) <= 1	):
			dx	=	abs(DetX[i] - DetX[i+1])
			dy	=	abs(DetY[i] - DetY[i+1])
			dist.append( 	((dx)**2	+	(dy)**2)**(0.5)	)
			match_1st.append(i)
			match_2nd.append(i+1)
	print("Same CZT N TICK events (within 0.040 s) = %d" %len(dist))
	#print(len(match_1st), len(match_2nd))

	dist_min	=	1.0
	dist_max	=	1.8
	ratio_min	=	1.0
	ratio_max	=	6.0
	
	## NON COMPTON CASE:    BUT WHY????????????????
	non_min_dist	=	2.0
	non_max_dist	=	100.0
	non_min_ratio	=	0.0	
	non_max_ratio	=	100.0
	
	Event_table_1st	=	cleaned_evt[match_1st]
	Event_table_2nd =	cleaned_evt[match_2nd]
	Distance	=	np.array(dist)
	
	## COMPTON EVENTS SORT OUT:
	adj_pix_stamp	=	[]
	remained_pix	=	[]
	for i in range(0, len(Event_table_1st)):
		if ((Distance[i] >= dist_min)	and	(Distance[i] <= dist_max)):
			adj_pix_stamp.append(i)
		else:
			remained_pix.append(i)
	print("Adjacent pixel events and at same time = %d" %len(adj_pix_stamp))

	adj_pix_1st	=	Event_table_1st[adj_pix_stamp]
	Time_1st	=	adj_pix_1st['Time']
	DetID_1st	=	adj_pix_1st['DetID']
	Energy_1st	=	adj_pix_1st['Energy']
	evt_dist	=	Distance[adj_pix_stamp]
	
	adj_pix_2nd	=	Event_table_2nd[adj_pix_stamp]
	Time_2nd	=	adj_pix_2nd['Time']
	DetID_2nd	=	adj_pix_2nd['DetID']
	Energy_2nd	=	adj_pix_2nd['Energy']
	##### We have 2 same time event pair one from Event_table_1st and other from Event_table_2nd.
	##### 2 compton event will be in same DetID.
	same_detid	=	[]
	for i in range(0, len(adj_pix_2nd)):
		if (	(DetID_1st[i]	==	DetID_2nd[i])	and (Energy_1st[i] > 0)	and	(Energy_2nd[i]	> 0)	):
			same_detid.append(i)
	print("Same det ID pairs = %d" % len(same_detid))
	adj_pix_table1	=	adj_pix_1st[same_detid]
	Time_1	=	adj_pix_table1['Time']
	DetID_1	=	adj_pix_table1['DetID']
	PixID_1	=	adj_pix_table1['pixID']
	Energy_1	=	adj_pix_table1['Energy']
	Tick_1	=	adj_pix_table1['cztnTick']
	DetX_1	=	adj_pix_table1['detX']
	DetY_1	=	adj_pix_table1['detY']
	dist_1_2	=	evt_dist[same_detid]

	adj_pix_table2	=	adj_pix_2nd[same_detid]
	Time_2	=	adj_pix_table2['Time']
	DetID_2	=	adj_pix_table2['DetID']
	PixID_2	=	adj_pix_table2['pixID']
	Energy_2	=	adj_pix_table2['Energy']
	Tick_2	=	adj_pix_table2['cztnTick']
	DetX_2	=	adj_pix_table2['detX']
	DetY_2	=	adj_pix_table2['detY']

	##### DECIDE WHICH PIX IS SCATTERED AND WHICH IS FOR ABSORBER IN EACH PAIR
	abs_1st_stamp	=	[]
	abs_2nd_stamp	=	[]
	for i in range(0, len(adj_pix_table1)):
		if (Energy_1[i]	>	Energy_2[i]):
			abs_1st_stamp.append(i)
		else:
			abs_2nd_stamp.append(i)		
	#print(len(abs_1st_stamp),len(abs_2nd_stamp))
	
	high_abs =	[]
	high_abs.extend(Energy_1[abs_1st_stamp])
	high_abs.extend(Energy_2[abs_2nd_stamp])
	
	low_scat	=	[]
	low_scat.extend(Energy_2[abs_1st_stamp])
	low_scat.extend(Energy_1[abs_2nd_stamp])
	
	high_abs	=	np.array(high_abs)
	low_scat	=	np.array(low_scat)
	total_energy	=	high_abs + low_scat
	ratio_h_by_l	=	high_abs/low_scat
	#print(len(ratio_h_by_l))
	
	Time	=	[]
	Time.extend(Time_1[abs_1st_stamp])
	Time.extend(Time_1[abs_2nd_stamp])
	Time	=	np.array(Time)
	#print(len(Time))
	evtdist	=	[]
	evtdist.extend(dist_1_2[abs_1st_stamp])
	evtdist.extend(dist_1_2[abs_2nd_stamp])
	evtdist	=	np.array(evtdist)
	#print(len(evtdist))
	h_detX	=	[]
	h_detX.extend(DetX_1[abs_1st_stamp])
	h_detX.extend(DetX_2[abs_2nd_stamp])
	h_detX	=	np.array(h_detX)
	h_detY	=	[]
	h_detY.extend(DetY_1[abs_1st_stamp])
	h_detY.extend(DetY_2[abs_2nd_stamp])
	h_detY	=	np.array(h_detY)
	l_detX	=	[]
	l_detX.extend(DetX_2[abs_1st_stamp])
	l_detX.extend(DetX_1[abs_2nd_stamp])
	l_detX	=	np.array(l_detX)
	l_detY	=	[]
	l_detY.extend(DetY_2[abs_1st_stamp])
	l_detY.extend(DetY_1[abs_2nd_stamp])	
	l_detY	=	np.array(l_detY)

	comp_eve	=	[]
	for i in range(0, len(Time)):
		if( 	( (ratio_h_by_l[i] >= ratio_min)	and	(ratio_h_by_l[i] < ratio_max) ) and	( (total_energy[i] >= args.emin) and (total_energy[i] <= args.emax) )	):
			comp_eve.append(i)
	print("Compton Events = %d" %len(comp_eve))
##=====================================================================================================================================
	Time	= 	Time[comp_eve]
	Q1_hist, bin_edges = np.histogram(Time, bins=tbins)
	#print(len(Q1_hist), max(Q1_hist))
	Q_max	=	Q_max + max(Q1_hist)
	## 12 Dec 2017	: LIVETIME correction is done: For 1 s lightcurve dividing the 1 s collected counts by fracexp of livetime fits.
##=====================================================================================================================================
	livetable	=	livefile[Q_n].data
	time	=	livetable['TIME']
	fracexp	=	livetable['FRACEXP']
	ltable	=	Table([time, fracexp], names=('Time', 'Fracexp'))
	#print(len(evtime))
	for i in range(0,len(time)):
		if (time[i] >= tmin):
			lstart = i
			break

	for j in range(0,len(time)):
		if (time[j] >= tmax):
			lstop = j
			break

	evtlive	=	ltable[lstart:lstop]
	ltime	=	evtlive['Time']
	lfrac	=	evtlive['Fracexp']
	Q1_hist	=	Q1_hist
	#Q1_hist	=	Q1_hist/lfrac
	#print ('Livetime correction done \n')
	Q_hist_sum = Q_hist_sum + Q1_hist
	#print(Q_hist_sum)
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0 
	#print(center_time[0],	center_time[-1])
Q_hist_sum_err	=	np.sqrt(Q_hist_sum)
#print(len(Q_hist_sum),len(Q_hist_sum_err))
## SELECT PRE BACKGROUND INTERVAL
print("\n SELECT PRE GRB BACKGROUND INTERVAL *********************")
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
plt.step(center_time, Q_hist_sum)
ax.set_xlim(tmin, tmax )
coords = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
pre_bkg_min=coords[0][0]
pre_bkg_max=coords[1][0]	
print("Pre-bkg start = %.2f , Pre-bkg stop = %.2f" %(pre_bkg_min, pre_bkg_max))

for i in range(0, len(center_time)):
	if(center_time[i]	>=	pre_bkg_min):
		stamp_prebkg_1	=  i
		break
for i in range(0, len(center_time)):
	if(center_time[i]	>=	pre_bkg_max):
		stamp_prebkg_2	=  i
		break

## SELECT GRB INTERVAL
print("\n SELECT THE GRB INTERVAL ********************************")
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
plt.step(center_time, Q_hist_sum)
ax.set_xlim(pre_bkg_max,tmax)
coords_grb = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
grb_start=coords_grb[0][0]
grb_stop=coords_grb[1][0]	
print("GRB start time = %.2f , GRB stop time = %.2f" %(grb_start, grb_stop))
for i in range(0, len(center_time)):
	if(center_time[i]	>=	grb_start):
		stamp_grb_1	=  i
		break
for i in range(0, len(center_time)):
	if(center_time[i]	>=	grb_stop):
		stamp_grb_2	=  i
		break

Q_hist_grb	=	Q_hist_sum[stamp_grb_1:stamp_grb_2]
#print(Q_hist_grb)
comp_sum	=	0.0
for i in range(0, len(Q_hist_grb)):
	comp_sum	=	comp_sum + Q_hist_grb[i]	
#print("Total Compton Events = %.2f " %(comp_sum))

## SELECT POST GRB INTERVAL:
print("\n SELECT THE POST GRB BACKGROUND INTERVAL ****************")
def onclick(event):
	global ix2, iy2
	ix2, iy2 = event.xdata, event.ydata
	global coords2
	coords2.append((ix1, iy1))
	if len(coords2) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(center_time,Q_hist_sum)
ax.set_xlim(grb_stop,tmax)
coords2 = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
post_bkg_min=coords2[0][0]
post_bkg_max=coords2[1][0]	
print("Post-bkg start = %.2f , Post-bkg stop = %.2f" %(post_bkg_min, post_bkg_max))
for i in range(0, len(center_time)):
	if(center_time[i]	>=	post_bkg_min):
		stamp_postbkg_1	=  i
		break
for i in range(0, len(center_time)):
	if(center_time[i]	>=	post_bkg_max):
		stamp_postbkg_2	=  i
		break

plt.errorbar(center_time, Q_hist_sum , yerr=Q_hist_sum_err,color='b',fmt="-*",ecolor='r',elinewidth=0.5,capsize=2)
plt.xlim((grb_start-100.0*args.tbin),(grb_stop+100.0*args.tbin))
plt.ylim(0,Q_max+14)
plt.xlabel('Time (In AstroSat time)',color='red')
plt.ylabel('Counts/s',color='red')
plt.title('Compton Event Lightcurve')
plt.axvline(center_time[stamp_grb_1], color='k', lw='1.10')
plt.axvline(center_time[stamp_grb_2], color='k',  lw='1.10')
plt.axvline(args.tmark, color='k', ls='--' ,lw='1.10')
plt.axvspan(grb_start, grb_stop, alpha=0.25,color='y')
plt.text((grb_start+2.0*args.tbin),(Q_max+0.1*Q_max), 'Total Compton Events = %.1f'% comp_sum, color='red',fontweight='bold')
if args.plottype == 'pdf':

	plotfile.savefig()
else:
	plt.show()	
	plt.savefig(args.outfile + "_compton_events." + args.plottype)
plt.show()	
if args.plottype == 'pdf':
    plotfile.close()


dblevtfile.close()
badfile.close()
livefile.close()

