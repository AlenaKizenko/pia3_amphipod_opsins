## We do not provide the initial file, as all the information is contained in the alignment (2.1_allopsins.aln)

## align
mafft --auto allopsins.faa >data/2.1_allopsins.aln ## L-INS-i
## trim
trimal -in data/2.1_allopsins.aln -out allopsins.trim.aln -automated1
## build the tree
iqtree -s allopsins.trim.aln -alrt 1000 -abayes -m LG+F+G4 # this model was the best option for almost all trees build by PIA3
