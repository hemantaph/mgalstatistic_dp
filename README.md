These are data products from the paper - More & More 2021 (see arxiv - TBD)

Below is the explanation for each directory.


1. AL - Advanced LIGO (3-detector - 1 year observing run)
2. Ap - Aplus (3-detector - 1 year observing run)
3. CE - Cosmic Explorer (single detector - 1 year long observing run)
4. ET - Einstein Telescope  (single detector - 1 year long observing run)

Each of the above directory contains -
a) lens catalog - ID, lens redshift, lens velocity dispersion (km/s), lens ellipticity, lens position angle (deg ), source redshift, source unlensed SNR, num_of_images (2-double, 4-quad) 
b) quadruply lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity 
c) doubly lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity 


5. O3 - LIGO-Virgo O3a observing run (3 detector - 6 months)

This directory contains (one extra column compared) -
a) lens catalog - ID, src_ID, lens redshift, lens velocity dispersion (km/s), lens ellipticity, lens position angle (deg ), source redshift, source unlensed SNR, num_of_images (2-double, 4-quad)  
b) quadruply lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity, src_ID 
c) doubly lensed image properties - ID, img_pos_x (arcsec), img_pos_y (arcsec), img_magnification, time_delay (days), parity, src_ID



Units and conventions:
* ID is the object identifier. It is the same between the lens and the image catalogs (double/quad) and helps to cross-match the same object between the two catalogs
* src_ID is an additional object identifier used only for O3 dir. Cross-match between the lens and image catalogs will require matching both ID and src_ID.
* Position angle : North is up - 0deg, East is left - 90deg
* Positions are relative to the lens galaxy at the center (0,0) in arcsec
* Magnification factors are absolute i.e. with respect to the source flux/intensity
* Time delays are relative to the image arriving first. The first image (type I) thus have a time delay of 0. 
* Parity: 1: positive parity corresponds to type I/ minimum  
  -1: negative parity corresponds to type II/saddle.
Note that type III/maximum are not formed for SIE lens models and these kind of images are not detected in this mock sample.

