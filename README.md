# HLS small area test

Test of HLS algorithm

!!! It is not official version of HLS !!!

Algorithm is in two parts. In first, atmospheric correction have to be processed (folder Atmospheric_correction). Then, the other parts (coregistration, resampling, BRDF-normalization, bandpass adjustment) are in scripts in folder HLS.

Input images are :
* Sentinel-2 : L1C
* Landsat 8 : Collection 1 Level 1
* Landsat 7 : Collection 1 Level 1
* Landsat 5 : Collection 1 Level 1

Input images for Landsat should be collection 2 but LaSRC for collection 2 is not available. It is the main problem of this algorithm.
For Landsat 7 and Landsat 5, it is only an attempt. Values for bandpass adjustment come from : D.P.Royet al. “Characterization of Landsat-7 to Landsat-8 reflective wavelength and norma-lized difference vegetation index continuity”

Results are in EPSG:32633.

In Test folder, there are scripts which enable to calculate differences between results.