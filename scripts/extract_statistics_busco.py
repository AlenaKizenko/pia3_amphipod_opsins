''' This is input example
# BUSCO version is: 3.0.2 
# The lineage dataset is: arthropoda_odb9 (Creation date: 2017-02-07, number of species: 60, number of BUSCOs: 1066)
# To reproduce this run: python /media/secondary/apps/busco/scripts/run_BUSCO.py -i /media/tertiary/opsins_data/Acanthogammarus_godlewskii.fsa_nt -o Acanthogammarus_godlewskii.fsa_nt_busco -l /media/tertiary/lineage/arthropoda_odb9/ -m transcriptome -c 1
#
# Summarized benchmarking in BUSCO notation for file /media/tertiary/opsins_data/Acanthogammarus_godlewskii.fsa_nt
# BUSCO was run in mode: transcriptome

	C:49.9%[S:40.9%,D:9.0%],F:20.5%,M:29.6%,n:1066

	532	Complete BUSCOs (C)
	436	Complete and single-copy BUSCOs (S)
	96	Complete and duplicated BUSCOs (D)
	219	Fragmented BUSCOs (F)
	315	Missing BUSCOs (M)
	1066	Total BUSCO groups searched'''
import os
#import csv
os.chdir('/home/alena/Documents/IB/project_opsins/busco/')
#with open('/home/alena/Documents/IB/project_opsins/stats.csv', 'w', newline='') as csvfile:
    #fieldnames = ['Species', 'Stats']
    #writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #writer.writeheader()
for filename in os.listdir('/home/alena/Documents/IB/project_opsins/busco/'):
    with open(filename) as file:
        gam = file.readlines()
        filename = filename[14:]
        #filename = filename[:-17]
        #writer.writerow({'Species': file_name, 'Stats': gam[7]})
        print(filename, gam[7]) # stdout was redirected to file "statistics.txt"
