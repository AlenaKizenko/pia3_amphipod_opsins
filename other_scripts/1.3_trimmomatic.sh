#!/bin/bash
lst=("SRR3532634" "SRR3532641" "SRR3532642") #listing reads` names
for i in ${!lst[@]};
do #performing Trimmomatic operations on listed reads
	name=${lst[$i]}
	echo $name
	textfastq1='_1.fastq.gz'
	textfastq2='_2.fastq.gz'
	textfastqout='.fastq'
	read1=$name$textfastq1
	read2=$name$textfastq2
	output=$name$textfastqout
	java -jar trimmomatic-0.36.jar PE ${PATH_TO_READS_DIR}/$read1 ${PATH_TO_READS_DIR}/$read2 -baseout ${PATH_TO_OUTPUT_DIR}/$output CROP:140 HEADCROP:20 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:36
done
