#!/bin/bash

  cat data/SRA_gammaridae.csv| while IFS=, read -r col1 col2
  do
	if [ -d "${PATH_TO_DIR}/${col2}" ] # if directory already exists
	then
		echo "Reads for $col2 are already downloaded"
	else
		mkdir ${PATH_TO_DIR}/${col2}
		echo "$col2 directory has been created"
		echo "$col1 reads start downloading"
    		fastq-dump --outdir ${PATH_TO_DIR}/$col2 --gzip --defline-seq '@$sn[_$rn]/$ri' --readids --split-3 $col1 # download reads via fastq-dump
		echo "$col1 reads have been downloaded in $col2 directory"
	fi
  done


