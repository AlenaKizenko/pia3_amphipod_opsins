#!/bin/bash

  cat /media/tertiary/Alena_Kizenko/SRA_gammaridae.csv| while IFS=, read -r col1 col2
  do
	if [ -d "/media/quartemary/rnaseq_reads/${col2}" ]
	then
		echo "Reads for $col2 are already downloaded"
	else
		mkdir /media/quartemary/rnaseq_reads/${col2}
		echo "$col2 directory has been created"
		echo "$col1 reads start downloading"
    		/media/secondary/apps/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump --outdir /media/quartemary/rnaseq_reads/$col2 --gzip --defline-seq '@$sn[_$rn]/$ri' --readids --split-3 $col1
		echo "$col1 reads have been downloaded in $col2 directory"
	fi
  done


