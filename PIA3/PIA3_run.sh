#!/bin/bash

  ls /media/tertiary/Alena_Kizenko/pia3_test | while read i;
  do
	python3 /media/tertiary/Alena_Kizenko/PIA3/PIA3.py -i /media/tertiary/Alena_Kizenko/pia3_test/${i} -o /media/tertiary/Alena_Kizenko/pia3_test/${i::-6}_pia3/ -db /media/tertiary/Alena_Kizenko/PIA3/classification_opsins_full_aa.fasta -s -all
  done
