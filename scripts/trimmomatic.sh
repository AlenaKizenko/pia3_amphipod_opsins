#!/bin/bash
lst=("SRR3532634" "SRR3532641" "SRR3532642")
for i in ${!lst[@]};
do
	name=${lst[$i]}
	echo $name
	textfastq1='_1.fastq.gz'
	textfastq2='_2.fastq.gz'
	textfastqout='.fastq'
	read1=$name$textfastq1
	read2=$name$textfastq2
	output=$name$textfastqout
	java -jar /media/secondary/apps/Trimmomatic-0.36/trimmomatic-0.36.jar PE /media/tertiary/Alena_Kizenko/hyalella_azteca_reads/$read1 /media/tertiary/Alena_Kizenko/hyalella_azteca_reads/$read2 -baseout /media/tertiary/Alena_Kizenko/hyalella_azteca_trim/$output CROP:140 HEADCROP:20 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:36
done
