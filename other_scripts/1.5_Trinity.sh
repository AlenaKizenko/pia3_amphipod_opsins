#!/bin/bash

  ls /media/quartemary/rnaseq_reads/Naumenko_reassemblies | while read i;
  do
	if [ -d "${PATH_TO_TRIMMED_READS_DIR}/trinity_$i.fasta" ]
	then
		echo "Transcriptome of $i is already assembled"
	else
		IFS=' ' read -r -a files <<< $(ls ${PATH_TO_TRIMMED_READS_DIR}/${i})
		read1=${files[0]::-9}_1P.fastq
		read2=${files[0]::-9}_2P.fastq
		mkdir ${PATH_TO_TRIMMED_READS_DIR}/trinity_$i
		Trinity --seqType fq --left ${PATH_TO_TRIMMED_READS_DIR}/${i}/${read1} --right ${PATH_TO_TRIMMED_READS_DIR}/${i}/${read2} --CPU 6 --max_memory 32G --full_cleanup --output ${PATH_TO_TRIMMED_READS_DIR}/trinity_$i
		echo ${i}
	fi
  done


