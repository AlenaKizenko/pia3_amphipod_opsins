#!/bin/bash

  ls /media/quartemary/rnaseq_reads/Naumenko_trimmed | while read i;
  do
	IFS=' ' read -r -a files <<< $(ls /media/quartemary/rnaseq_reads/Naumenko_trimmed/${i})
	read1=${files[0]::-9}_1P.fastq
	read2=${files[0]::-9}_2P.fastq
	mkdir /media/tertiary/Alena_Kizenko/rnaspades/${i}_rnaspades
	/media/secondary/apps/SPAdes-3.13.1-Linux/bin/rnaspades.py -t 8 -1 /media/quartemary/rnaseq_reads/Naumenko_trimmed/${i}/${read1} -2 /media/quartemary/rnaseq_reads/Naumenko_trimmed/${i}/${read2} -o /media/tertiary/Alena_Kizenko/rnaspades/${i}_rnaspades
	echo ${i}
  done
