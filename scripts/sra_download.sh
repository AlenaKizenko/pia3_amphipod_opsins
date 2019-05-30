#!/bin/bash

lst=("SRR3532634" "SRR3532641" "SRR3532642")
for i in ${!lst[@]};
do
	/media/secondary/apps/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump --outdir /media/tertiary/Alena_Kizenko/hyalella_azteca_reads --gzip --defline-seq '@$sn[_$rn]/$ri' --readids --split-files ${lst[$i]}
done
