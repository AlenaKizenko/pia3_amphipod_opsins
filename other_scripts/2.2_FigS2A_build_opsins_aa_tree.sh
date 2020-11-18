mafft --auto allopsins.faa >allopsins.aln ##L-INS-i
$appdir/trimal/source/trimal -in allopsins.aln -out allopsins.trim.aln -automated1
$appdir/iqtree-1.6.10-Linux/bin/iqtree -s allopsins.trim.aln -alrt 1000 -abayes -m LG+F+G4

