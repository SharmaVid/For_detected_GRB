#!/bin/bash

# for_detected_grb_noise_pol.sh
# Priyanka & Vidushi, Jan, 2018
 
# Version 1.0------------------------------------------
### It creates three folders for detected GRBs in directory: special/GRBname
### (1) latest_modeM0_up:  It generates quad_clean.evt, given time bin livetime.fits and cztbindata lc. Here calculate_T90.py code is run for calculation of several parameters like T90, Total count , peak count.
### (2) polarization_analysis: It creates directory for polarization analysis for threshold 200 and 4. The compton events counts and lightcurve is generated here. 
### (3) noise_reduced: 


echo "####################################  Parsing the argument for processing Gamma Ray Burst and Finding the orbit containing Gamma Ray Burst ##############################################"
GRBNAME=$1
#DATE=$2
#TIME=$3
Astrosat_seconds=$2
RA=$3
DEC=$4
Tbin=$5
#orbit_num=`Finding_orbit.sh $DATE $TIME`
#grblogs=/data2/czti/special/grblogs/$GRBNAME"_after_GRB_detected.log"
if [ -z "$5" ]
  then
  	Tbin=1  
fi



echo "####################################  Converting UTC to AstroSat time ##############################################"
#Time=`echo $DATE $TIME`
#Output=$( Astrosat_Time.py -iso "$Time")
#Astrosat_seconds=`echo $Output|cut -d ',' -f2`
echo "#########################################  Defining Defalut path #######################################"
specail=/data2/czti/special
noise_reduced=/home/cztipoc/czti/users/ajayr/noise_codes/config
noise_reduced1=/home/cztipoc/czti/users/ajayr/noise_codes/bin
caldbbadpix=/home/cztipoc/CALDB/data/as1/czti/bcf
mkf=$(find $specail/$GRBNAME/ -name "*.mkf")
echo "MKF file: $mkf"
mkf1=$(basename $mkf)
orbit_num=`echo $mkf1 |cut -d '_' -f4|cut -d 'c' -f1`
echo -e ' \t ' 
echo "######################################  Start Processing $GRBNAME  #######################################"
echo "######################################  Start Processing $GRBNAME  #######################################" #>>$grblogs
echo "#####Processing Time:: $(date "+%Y-%m-%d %H:%M:%S")" #>>$grblogs

echo -e ' \t' 
echo "#$GRBNAME corresponds to orbit:$orbit_num#"
echo "$GRBNAME corresponds to orbit:$orbit_num" #>>$grblogs
echo -e ' \t'

echo "Input Provided:: 
 GRBNAME:          $GRBNAME 
 Orbit Number:     $orbit_num 
 " #>>$grblogs 



echo "Input Provided:: 
 GRBNAME:          $GRBNAME 
 Orbit Number:     $orbit_num 
 "
 


tmark=$Astrosat_seconds
tmin=`echo "$Astrosat_seconds-500"|bc -l`
tmax=`echo "$Astrosat_seconds+500"|bc -l`






echo -e ' \t'
echo "#####################################  Creating directory Structure for $GRBNAME #########################"
echo "#####Creating directory Structure for $GRBNAME:: " #>>$grblogs
echo -e ' \t'
echo "`echo -e ' \t'`" #>>$grblogs
     
  
     mkdir $specail/$GRBNAME/noise_reduced
     echo "mkdir $specail/$GRBNAME/noise_reduced" #>>$grblogs
     mkdir $specail/$GRBNAME/pol_analysis
     echo "mkdir $specail/$GRBNAME/pol_analysis" #>>$grblogs
     mkdir $specail/$GRBNAME/latest_modeM0_up
     echo "mkdir $specail/$GRBNAME/latest_modeM0_up" #>>$grblogs
     
    #chmod 755 $specail/$GRBNAME/latest_modeM0_up	
  
echo "##################################  Defining new path and variable name ########################################################"
echo "#####Defining new path and variable name::"  #>>$grblogs
echo "`echo -e ' \t'`"  #>>$grblogs


   modeM0=$specail/$GRBNAME/modeM0
   latest_modeM0=$specail/$GRBNAME/latest_modeM0_up
   NR_products=$specail/$GRBNAME/noise_reduced
   Pol_products=$specail/$GRBNAME/pol_analysis
   echo "modeM0=$specail/$GRBNAME/modeM0" #>>$grblogs
   echo "Pol_products=$specail/$GRBNAME/pol_analysis" #>>$grblogs



echo "#############################################  Finding Products for $GRBNAME   #################################################"
   
   bcevt=$(find $modeM0/ -name "*level2_bc.evt")
   echo "Bunch Clean event file: $bcevt"
   bunch_clean_livetime=$(find  $modeM0/ -name "*level2_bc_livetime.fits")
   echo "Bunch Clean livetime file: $ bunch_clean_livetime"
    evt=$(find $modeM0/ -name "*level2_bc.evt")
   echo "Bunch Clean event file: $evt"

   bunnchevt=$(find $modeM0/ -name "*cztM0_level2_bunch.fits")
   echo $bunnchevt

    cp $bcevt $latest_modeM0/
    cp $bunch_clean_livetime $latest_modeM0/
    cp $bcevt $Pol_products/
    cp $bunch_clean_livetime $Pol_products/
    cp $evt $NR_products/
    cp $bunnchevt $NR_products/



	quad_bc_ds_evt=`echo $bcevt|sed 's/bc.evt/quad_bc_ds.evt/'`
	quad_bc_ds_evt=$latest_modeM0/$(basename $quad_bc_ds_evt)


 echo $quad_bc_ds_evt

echo "####################################  Genrating Datasel file for processing for  Calculating T90  ###########################"
echo "#####Genrating Datasel file for processing::" #>>$grblogs
echo "`echo -e ' \t'`" #>>$grblogs


         cztdatasel infile=$bcevt gtifile=$bcevt gtitype="QUAD" outfile=$quad_bc_ds_evt clobber="y" history="y"
         echo "cztdatasel infile=$bcevt gtifile=$bcevt gtitype="QUAD" outfile=$quad_bc_ds_evt clobber="y" history="y"" #>>$grblogs
 
 
 	echo $quad_bc_ds_evt
	quad_pc_evt=`echo $quad_bc_ds_evt |sed 's/quad_bc_ds.evt/quad_bc_ds_pc.evt/'`
	quad_pc_evt=$latest_modeM0/$(basename $quad_pc_evt)
	echo $quad_pc_evt

	#quad_pc_double_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_bc_ds_pc.dblevt/'`
	#quad_pc_double_evt=$latest_modeM0/$(basename $quad_pc_double_evt)
	#echo $quad_pc_double_evt

	quad_livetime=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_livetime.fits/'`
	quad_livetime=$latest_modeM0/$(basename $quad_livetime)
	echo $quad_livetime

	quad_badpix=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_badpix.fits/'`
	quad_badpix=$latest_modeM0/$(basename $quad_badpix)
	echo $quad_badpix

echo -e ' \t'
echo "####################################  Genrating pixclean evt file for processing  Calculating T90 ###########################"
echo "#####Genrating pixclean evt file for processing::" #>>$grblogs
echo -e ' \t'
echo "`echo -e ' \t'`" #>>$grblogs


            echo "cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix  par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=1000 par_pix_count_thresh=100" #>>$grblogs


            cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix  par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=1000 par_pix_count_thresh=100

echo "###################################   Genrating common clean evt file  Calculating T90 #######################################"
echo "#####Genrating common clean evt file::" #>>$grblogs
echo "`echo -e ' \t'`" #>>$grblogs
  
           quad_clean_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_clean.evt/'`
            echo "cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"" #>>$grblogs 
	   quad_clean_evt=$latest_modeM0/$(basename $quad_clean_evt)
	   cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"
	
	
	   #quad_clean_dblevt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_clean.dblevt/'`
	   #echo " cztevtclean infile=$quad_pc_double_evt outfile=$quad_clean_dblevt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"" #>>$grblogs  
	   #quad_clean_dblevt=$latest_modeM0/$(basename $quad_clean_dblevt)
          # cztevtclean infile=$quad_pc_double_evt outfile=$quad_clean_dblevt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"

echo "###################################   Genrating cztbindata lc file #######################################"
echo "#####Genrating common clean evt file::" #>>$grblogs
echo "`echo -e ' \t'`" #>>$grblogs

outfile=`echo $quad_badpix |sed 's/level2_quad_badpix.fits/level2_grb/'`
outfile=$latest_modeM0/$(basename $outfile)
outevtfile=`echo $quad_badpix|sed 's/level2_quad_badpix.fits/level2_weight_grb.evt/'`
outevtfile=$latest_modeM0/$(basename $outevtfile)

cztbindata inevtfile=$quad_clean_evt  mkffile=$mkf badpixfile=$quad_badpix livetimefile=$quad_livetime outfile=$outfile outevtfile=$outevtfile maskWeight="no" rasrc=$RA decsrc=$DEC badpixThreshold="0" outputtype="lc"

calculate_T90="AS1CZT_"$GRBNAME"_"
echo "################################### Calculating T90 #######################################"
echo "#####Calculating T90::" #>>$grblogs
echo "`echo -e ' \t'`" #>>$grblogs
  
  python /data2/czti/special/grb_detection_script/calculate_T90_V1.py $quad_clean_evt $caldbbadpix/AS1cztbadpix20160908v01.fits $quad_livetime $latest_modeM0/$calculate_T90 $GRBNAME  --tmark $Astrosat_seconds 
echo "Do you want to Run Noise Reduced pipeline????"
read Run
if [ "$Run" == "yes" ] || [ "$Run" == "y" ] || [ "$Run" == "Yes" ] ||  [ "$Run" == "Y" ] ||  [ "$Run" == "YES" ]   ; then
	echo "####################################  Genrating NOISE Reduced file for processing ###########################"
	echo "##### Genrating Noice Reduced file for processing::" #>>$grblogs
	echo "`echo -e ' \t'`" #>>$grblogs

	bash $noise_reduced1/clean_data.sh $mkf $noise_reduced/saaThreshold $NR_products/$(basename $evt) $noise_reduced/AS1cztbadpix20160908v01.fits $NR_products/$(basename $bunnchevt) $noise_reduced/noiseReductionThreshold $noise_reduced/AS1cztlld20160517v01.fits $bunnchevt 
	mv *cztM0_level2_*_* $NR_products
        mv *cztM0_level2_*.evt $NR_products
        mv *cztM0_level2_ds.log $NR_products

        #echo "################################### Calculating T90 for noise reduced #######################################"
        #echo "#####Calculating T90::" #>>$grblogs
        #echo "`echo -e ' \t'`" #>>$grblogs
           #python calculate_T90.py $quad_clean_evt $caldbbadpix/AS1cztbadpix20160908v01.fits $quad_livetime $GRBNAME --tmark $Astrosat_seconds 

fi
       echo "####################################  Genrating Datasel file for processing for Generating polarisation analysis  ###########################"
       echo "#####Genrating Datasel file for processing::" #>>$grblogs
       echo "`echo -e ' \t'`" #>>$grblogs

        quad_bc_ds_evt=`echo $bcevt|sed 's/bc.evt/quad_bc_ds.evt/'`
	quad_bc_ds_evt=$Pol_products/$(basename $quad_bc_ds_evt)
        echo $quad_bc_ds_evt

        cztdatasel infile=$bcevt gtifile=$bcevt gtitype="QUAD" outfile=$quad_bc_ds_evt clobber="y" history="y"
        echo "cztdatasel infile=$bcevt gtifile=$bcevt gtitype="QUAD" outfile=$quad_bc_ds_evt clobber="y" history="y"" #>>$grblogs
 
 
 	echo $quad_bc_ds_evt
	quad_pc_evt=`echo $quad_bc_ds_evt |sed 's/quad_bc_ds.evt/quad_bc_ds_pc.evt/'`
	quad_pc_evt=$Pol_products/$(basename $quad_pc_evt)
	echo $quad_pc_evt

	quad_pc_double_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_bc_ds_pc.dblevt/'`
	quad_pc_double_evt=$Pol_products/$(basename $quad_pc_double_evt)
	echo $quad_pc_double_evt

	quad_livetime=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_livetime.fits/'`
	quad_livetime=$Pol_products/$(basename $quad_livetime)
	echo $quad_livetime

	quad_badpix=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_badpix.fits/'`
	quad_badpix=$Pol_products/$(basename $quad_badpix)
	echo $quad_badpix

	quad_pc_evt1=`echo $quad_bc_ds_evt |sed 's/quad_bc_ds.evt/quad_bc_ds_pc.evt/'`
	quad_pc_evt1=$Pol_products/$(basename $quad_pc_evt)
	echo $quad_pc_evt1

	quad_pc_double_evt1=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_bc_ds_pc.dblevt/'`
	quad_pc_double_evt1=$Pol_products/$(basename $quad_pc_double_evt)
	echo $quad_pc_double_evt1

	quad_livetime1=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_livetime.fits/'`
	quad_livetime1=$Pol_products/$(basename $quad_livetime)
	echo $quad_livetime1

	quad_badpix1=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_badpix.fits/'`
	quad_badpix1=$Pol_products/$(basename $quad_badpix)
	echo $quad_badpix1

       echo -e ' \t'
       echo "####################################  Genrating pixclean evt file for processing Generating polarisation analysis with Threshold values 200/4 & user given ratio  ###########################"
       echo "#####Genrating pixclean evt file for processing::" #>>$grblogs
       echo -e ' \t'
       echo "`echo -e ' \t'`" #>>$grblogs


            echo "cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix par_outfile2=$quad_pc_double_evt par_writedblevt="y" par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=200 par_pix_count_thresh=4" #>>$grblogs


            cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix par_outfile2=$quad_pc_double_evt par_writedblevt="y" par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=200 par_pix_count_thresh=4
   
            echo "Enter par_det_count_thresh and par_pix_count_thresh::"
             read det_th pix_th
	
	    echo "cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt1 par_outlivetimefile=$quad_livetime1  par_badpixfile=$quad_badpix1 par_outfile2=$quad_pc_double_evt1 par_writedblevt="y" par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=$det_th par_pix_count_thresh=$pix_th" #>>$grblogs


            cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$bunch_clean_livetime par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix par_outfile2=$quad_pc_double_evt par_writedblevt="y" par_nsigma=5 par_det_tbinsize=$Tbin par_pix_tbinsize=$Tbin par_det_count_thresh=$det_th par_pix_count_thresh=$pix_th


        echo "###################################   Genrating common clean evt file  Generating polarisation analysis#######################################"
        echo "#####Genrating common clean evt file::" #>>$grblogs
        echo "`echo -e ' \t'`" #>>$grblogs
  
           quad_clean_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_clean.evt/'`
            echo "cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"" #>>$grblogs 
	   quad_clean_evt=$Pol_products/$(basename $quad_clean_evt)
	   cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"
	
	
	   quad_clean_dblevt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_clean.dblevt/'`
	   echo " cztevtclean infile=$quad_pc_double_evt outfile=$quad_clean_dblevt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"" #>>$grblogs  
	   quad_clean_dblevt=$Pol_products/$(basename $quad_clean_dblevt)
           cztevtclean infile=$quad_pc_double_evt outfile=$quad_clean_dblevt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"

          Pol_compton_events="AS1CZT_"$GRBNAME"_"

        echo "################################### Generating polarisation analysis data #######################################"
        echo "#####creating polarisation analysis data::" #>>$grblogs 
        echo "`echo -e ' \t'`" #>>$grblogs

         python /data2/czti/special/grb_detection_script/compton_events.py $quad_clean_dblevt $caldbbadpix/AS1cztbadpix20160908v01.fits  $quad_livetime $Pol_products/$Pol_compton_events --tmark $Astrosat_seconds 

        echo "###################################   Deleting Uncessary outputs #######################################"
        echo "#####Deleting Uncessary outputs" #>>$grblogs
        echo "`echo -e ' \t'`" #>>$grblogs


           delete=$(find $latest_modeM0/ -name  *bc_ds_pc.evt -o -name *bc_livetime.fits -o -name *quad_bc_ds.evt -o -name *quad_badpix.fits -o -name *bc.evt )
           echo "The following file is being deleted: $delete"
           echo "The following file is being deleted: $delete" #>>$grblogs
           rm  -rf $delete

          delete=$(find $Pol_products/ -name *quad_bc_ds_pc.dblevt -o -name *bc_ds_pc.evt -o -name *bc_livetime.fits -o -name *quad_bc_ds.evt -o -name *quad_badpix.fits -o -name *bc.evt )
           echo "The following file is being deleted: $delete"
           echo "The following file is being deleted: $delete" #>>$grblogs
           rm -rf $delete



	

