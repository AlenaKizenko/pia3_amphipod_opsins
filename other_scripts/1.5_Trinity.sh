#!/bin/bash

  ls /media/quartemary/rnaseq_reads/Naumenko_reassemblies | while read i;
  do
	if [ -d "/media/tertiary/Alena_Kizenko/Naumenko_reassemblies/trinity_$i.fasta" ]
	then
		echo "Transcriptome of $i is already assembled"
	else
		IFS=' ' read -r -a files <<< $(ls /media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i})
		read1=${files[0]::-9}_1P.fastq
		read2=${files[0]::-9}_2P.fastq
		mkdir /media/tertiary/Alena_Kizenko/Naumenko_reassemblies/trinity_$i
		/media/secondary/apps/trinityrnaseq-2.8.5/Trinity --seqType fq --left /media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i}/${read1} --right /media/quartemary/rnaseq_reads/Naumenko_reassemblies/${i}/${read2} --CPU 6 --max_memory 32G --full_cleanup --output /media/tertiary/Alena_Kizenko/Naumenko_reassemblies/trinity_$i
		echo ${i}
	fi
  done


