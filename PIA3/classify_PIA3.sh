#!/bin/bash

  ls /media/tertiary/Alena_Kizenko/merged_PIA3/ | while read i;
  do
	cd /media/tertiary/Alena_Kizenko/merged_PIA3/${i}	
	
	path=$PWD

	file=${i::-5}

	python3 /media/tertiary/Alena_Kizenko/PIA3/class_opsins.py query_class_align.fasta.contree $file $path
	
	cp ${path}/${file}_opsins_class.fasta /media/tertiary/Alena_Kizenko/opsins_merged_aa_PIA3/${file}_opsins_class_PIA3.fasta

	#cp ${path}/${file}_opsins_class.csv /media/tertiary/Alena_Kizenko/opsins_merged_aa_PIA3_csv/${file}_opsins_class_PIA3.csv
	
	echo ${file}

	cd ..

  done

