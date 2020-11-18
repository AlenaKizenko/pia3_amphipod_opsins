#!/bin/bash

  cat amphipoda_TSA.tsv| while IFS=$'\t' read -r col1 col2
	mkdir /media/quartemary/TSA/${col2}_${col1}
	echo "$col2 $col1 started downloading"
    	/media/secondary/apps/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump --outdir /media/quartemary/TSA/${col2}_${col1} -F --fasta $col1
	echo "$col2 $col1 ready"
  done

## downloads all amphipod assemblies from TSA (99 as of 16 Oct 2019)

