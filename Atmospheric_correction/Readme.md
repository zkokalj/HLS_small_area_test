# Atmospheric correction

For Landsat 8 and Sentinel-2, LaSRC 2.0.1 is used. For Landsat 5 and Landsat 7, LEDAPS 3.4.0 is used.
Be careful : in HLS 1.5, version of LaSRC is 3.5.5. But this version is not currently available (2021/07/27). The main inconvenience is that only in version 2.0.1 collection 1 level 1 images are accepted, not collection 2 level 1.
Algorithm comes from here : https://github.com/marujore/LaSRC-LEDAPS-Fmask
run.sh enables to run this algorithm on all images which are in the same folder.

# Structure of files
* run.sh
* LaSRC : contains main scripts
    * Dockerfile : install on the computer LaSRC and LEDAPS codes
    * Fmask_4_3_Linux.install : process Fmask algorithm
    * run_lasrc_ledaps_fmask.sh : process LaSRC and LEDAPS algorithm
* LaSRC_auxiliary : auxiliary data for Landsat 8 and Sentinel-2 : https ://edclpdsftp.cr.usgs.gov/downloads/auxiliaries/lasrc_auxiliary/L8/
* Ledaps_auxiliary : auxiliary data for Landsat 5 and Landsat 7 : https ://edclpdsftp.cr.usgs.gov/downloads/auxiliaries/ledaps_auxiliary/L17/
* Input_images : 
    * S2....SAFE (for Sentinel-2)
    * LC08... (for Landsat 8)
    * LE07... (for Landsat 7)
    * LT05... (for Landsat 5)
* Results : result of algorithm