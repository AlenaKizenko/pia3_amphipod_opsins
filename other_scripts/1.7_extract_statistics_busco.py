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


os.chdir('PATH_TO_BUSCO_STATISTICS')
for filename in os.listdir('PATH_TO_BUSCO_STATISTICS'):
    with open(filename) as file:
        gam = file.readlines()
        filename = filename[14:] # get amphipode name from file name
        print(filename, gam[7]) # get line with statistics
	# stdout was redirected to file "statistics.txt"
