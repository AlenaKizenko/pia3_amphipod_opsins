#!/bin/bash

  ls $PATH_TO_TRIMMED_READS_DIR | while read i;
  do
	IFS=' ' read -r -a files <<< $(ls ${PATH_TO_TRIMMED_READS_DIR}/${i})
	read1=${files[0]::-9}_1P.fastq
	read2=${files[0]::-9}_2P.fastq
	mkdir ${PATH_TO_TRIMMED_READS_DIR}/${i}_rnaspades
	rnaspades.py -t 8 -1 ${PATH_TO_TRIMMED_READS_DIR}/${i}/${read1} -2 ${PATH_TO_TRIMMED_READS_DIR}/${i}/${read2} -o ${PATH_TO_TRIMMED_READS_DIR}/${i}_rnaspades
	echo ${i}
  done
