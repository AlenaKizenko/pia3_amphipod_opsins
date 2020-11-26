#!/bin/bash

  cat 1.2_amphipoda_TSA.tsv | while IFS=$'\t' read -r col1 col2
	mkdir ${PATH_TO_DIR}/${col2}_${col1}
	echo "$col2 $col1 started downloading"
    	fastq-dump --outdir ${PATH_TO_DIR}/${col2}_${col1} -F --fasta $col1
	echo "$col2 $col1 ready"
  done

## downloads all amphipod assemblies from TSA (99 as of 16 Oct 2019)

