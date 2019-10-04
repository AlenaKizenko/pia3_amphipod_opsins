#!/bin/bash

  ls /media/quartemary/rnaseq_reads/ | while read i;
  do
	if [ -d "/media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i}" ]
	then
		echo "Reads for ${i} are already trimmed"
	else
		IFS=' ' read -r -a files <<< $(ls /media/quartemary/rnaseq_reads/${i})
		fq=".fastq"
		output=${files[0]::-11}$fq
		mkdir /media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i}
		java -jar /media/secondary/apps/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 12 /media/quartemary/rnaseq_reads/${i}/${files[0]} /media/quartemary/rnaseq_reads/${i}/${files[1]} -baseout /media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i}/${output} ILLUMINACLIP:/media/secondary/apps/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 MINLEN:50 AVGQUAL:20
		echo $i
	fi
  done


