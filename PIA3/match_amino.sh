#!/bin/bash

  ls /media/tertiary/Alena_Kizenko/merged_PIA3 | while read i;
  do
	cd /media/tertiary/Alena_Kizenko/merged_PIA3/$i	
	
	path=$PWD

	file=${i::-5}

	python3 /media/tertiary/Alena_Kizenko/PIA3/match_amino_nucl.py $path $file

	cp ${path}/${file}_opsins_nucl.fasta /media/tertiary/Alena_Kizenko/opsins_merged_nucl_PIA3/${file}_opsins_nucl.fasta
	
	echo ${file}

	cd ..

  done

