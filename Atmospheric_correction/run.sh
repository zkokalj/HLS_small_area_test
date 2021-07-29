#!/bin/bash

cheminSentinel=`pwd`/Input_images:/mnt/input-dir:rw
cheminL8=`pwd`/LaSRC_auxiliary:/mnt/lasrc-aux:ro
cheminL57=`pwd`/Ledaps_auxiliary:/mnt/ledaps-aux:ro
cheminResults=`pwd`/Results:/mnt/output-dir:rw

cd Input_images

liste=`ls`


for fichier in $liste
do 
	echo $fichier
	if [[ $fichier == "LT05"* ]] || [[ $fichier == "LE07"* ]]; then
	docker run --rm -v $cheminSentinel -v $cheminResults -v $cheminL57 -t lasrc_ledaps_fmask $fichier
	else
	docker run --rm -v $cheminSentinel -v $cheminResults -v $cheminL8 -t lasrc_ledaps_fmask $fichier
	fi
done





