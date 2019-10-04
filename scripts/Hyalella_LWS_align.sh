#!/bin/bash

  ls /media/quartemary/rnaseq_reads/ | while read i;
  do
	IFS=' ' read -r -a files <<< $(ls /media/quartemary/rnaseq_reads/${i})
	/media/secondary/apps/bowtie2-2.3.2/bowtie2 -x /media/tertiary/Alena_Kizenko/Hyalella_azteca_LWS/index_HLWS_opsins -p 12 -1 /media/quartemary/rnaseq_reads/${i}/${files[0]} -2 /media/quartemary/rnaseq_reads/${i}/${files[1]} -S /media/tertiary/Alena_Kizenko/Hyalella_azteca_LWS/${i}.sam &> /media/tertiary/Alena_Kizenko/Hyalella_azteca_LWS/${i}.txt
	rm /media/tertiary/Alena_Kizenko/Hyalella_azteca_LWS/${i}.sam
	echo ${i} ${files[0]} ${files[1]}
  done


