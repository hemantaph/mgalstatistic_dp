These are data products from the paper - More & More 2021 (see arxiv - TBD)

A) Summary:

1. mkplots\*.py -  will generate figures from the paper except for the O3_SIS model using the catalogs and scripts listed below.

2. snr\*.py -
If SNRthresh is set to 8 (default), it allows one to choose only those pairs from a given pair combination (e.g. 31 ) where the fainter lensed image (e.g. image 3) is above networkSNR of 8.
If SNRthresh is set to 0 in snr\*.py, then all those quads where the 2nd brightest lensed image will have networkSNR > 8 are used in the calculation.

3. a) unlensedpairs - comprise of the distributions for the unlensed pair population
   
   b) gwevents_O3 - comprise of info on event pairs taken from Table 3 of Abbott et al. 2021 (o3a lensing paper)
   
   c) out_fincat\_\* - comprise of the distributions for the lensed pair population for various detectors (see further details below). 




B) Below is the explanation for each directory for the lens population

1. out_fincat_AL - Advanced LIGO (3-detector - 1 year observing run)
2. out_fincat_Ap - Aplus (3-detector - 1 year observing run)
3. out_fincat_CE - Cosmic Explorer (single detector - 1 year long observing run)
4. out_fincat_ET - Einstein Telescope  (single detector - 1 year long observing run)

   Each of the above directory contains -

   a) lens catalog - ID, lens redshift, lens velocity dispersion (km/s), lens ellipticity, lens position angle (deg ), source redshift, source unlensed SNR, num_of_images (2-double, 4-quad) 

   b) quadruply lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity 
   
   c) doubly lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity 


5. out_fincat_O3 - LIGO-Virgo O3a observing run (3 detector - 6 months)

   This directory contains (one extra column compared) -
   
   a) lens catalog - ID, src_ID, lens redshift, lens velocity dispersion (km/s), lens ellipticity, lens position angle (deg ), source redshift, source unlensed SNR, num_of_images (2-double, 4-quad)  
   
   b) quadruply lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity, src_ID 
   
   c) doubly lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity, src_ID



Definitions, units and conventions:
* ID is the object identifier. It is the same between the lens and the image catalogs (double/quad) and helps to cross-match the same object between the two catalogs
* src_ID is an additional object identifier used only for O3 dir. Cross-match between the lens and image catalogs will require matching both ID and src_ID.
* Position angle is measured counter-clockwise in deg where North is up - 0deg, East is left - 90deg
* Positions are in arcsec relative to the lens galaxy at the center (0,0)
* Magnification factors are absolute meaning they are with respect to the source flux/intensity. In other words, unlensed SNR x magnification factor gives the lensed image SNR given the sensitivity of a detector
* Time delays are given in days relative to the image arriving first. The first image (type I) thus have a time delay of 0. 
* Parity: 1: positive parity corresponds to type I/ minimum  
  -1: negative parity corresponds to type II/saddle
Note that type III/maximum are not formed for SIE lens models and these kind of images are not detected in this mock sample.

