#!/bin/bash

  ls /media/tertiary/Alena_Kizenko/rest_pia3 | while read i;
  do
	python3 /media/tertiary/Alena_Kizenko/PIA3/PIA3.py -i /media/tertiary/Alena_Kizenko/rest_pia3/${i} -o /media/tertiary/Alena_Kizenko/merged_PIA3/${i::-6}_pia3/ -db /media/tertiary/Alena_Kizenko/PIA3/classification_opsins_full_aa.fasta -s -cds
  done
