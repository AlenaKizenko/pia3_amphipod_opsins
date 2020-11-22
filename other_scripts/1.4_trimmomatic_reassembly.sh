#!/bin/bash

  ls $PATH_TO_READS_DIR | while read i;
  do
	if [ -d "${PATH_TO_TRIMMED_READS_DIR}/${i}" ]
	then
		echo "Reads for ${i} are already trimmed"
	else
		IFS=' ' read -r -a files <<< $(ls ${PATH_TO_READS_DIR}/${i})
		fq=".fastq"
		output=${files[0]::-11}$fq
		mkdir ${PATH_TO_TRIMMED_READS_DIR}/${i}
		java -jar trimmomatic-0.36.jar PE -threads 12 ${PATH_TO_READS_DIR}/${i}/${files[0]} ${PATH_TO_READS_DIR}/${i}/${files[1]} -baseout ${PATH_TO_TRIMMED_READS_DIR}/${i}/${output} ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 MINLEN:50 AVGQUAL:20
		echo $i
	fi
  done


