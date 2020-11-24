## Preamble: 
## All assemblies have too few common one-copy orthologs, and making a presence/absence template is too complicated.
## Also, we should better avoid bias in taxonomic sampling.
## Thus, we would need to choose an adequate number of assemblies to proceed.
## Criteria: 
## * only one species per genus (excluding the formal genus Gammarus; there, only one per species);
## * only assemblies with > 50% BUSCO completeness (original TSA assembly);
## * only assemblies with at least one opsin found.

## Get all necessary assemblies.
## Reproducibility: the assemblies are available from Dryad.

#export ${DIR_WITH_ASSEMBLIES}=/enter/your/path/here

cp ${DIR_WITH_ASSEMBLIES}/Acanthogammarus_godlewskii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Asprogammarus_rhodophthalmus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Boeckaxelia_carpenterii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Brandtia_latissima_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Caprella_sp_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Carinurus_bicarinatus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Dorogostaiskia_parasitica_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Echinogammarus_berilloni_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Echiuropus_macronychus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Eulimnogammarus_cyaneus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_chevreuxi_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_fossarum_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_lacustris_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_minus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_pisinnus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_pulex_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gammarus_wautieri_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Gmelinoides_fasciatus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Heterogammarus_sophianosii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Hirondellea_gigas_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Homalogammarus_brandtii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Hyalellopsis_costata_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Linevichella_vortex_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Macropereiopus_parvus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Marinogammarus_marinus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Micruropus_wahlii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Odontogammarus_calcaratus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Ommatogammarus_albinus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Oxyacanthus_curtus_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Pachyschesis_branchialis_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Palicarinus_puzyllii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Pallasea_cancelloides_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Parapallasea_wosnessenskii_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Parhyale_hawaiensis_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Poekilogammarus_pictoides_rnaspades.fasta .
cp ${DIR_WITH_ASSEMBLIES}/Talitrus_saltator_rnaspades.fasta .

## Predict proteins with TransDecoder
for assembly in *_rnaspades.fasta; do TransDecoder.LongOrfs -t $assembly; TransDecoder.Predict -t $assembly --single_best_only; done

## Cluster very similar proteins ## .95 is the best option according to the analysis of G. minus opsins
for opsinlist in *pep; do cd-hit -i $opsinlist -c .95 -o $opsinlist.repr.faa; done
mkdir cdhit_aa
mv *faa cdhit_aa/
mv *clstr cdhit_aa/

## Build orthologous groups
proteinortho6.pl  *faa -cpus=6 -project=36amph


			```{r}
			options(stringsAsFactors = F)
		
			## read the proteinortho table
			tbl <- read.delim("36amph.proteinortho", na.strings = "*")
					
			## one-to-one orthologs for the tree!
			## present in each species
			presentinall <- tbl[complete.cases(tbl), ]
			## present in each species only once ## in this case each sample
			onetoone <- presentinall[presentinall$X..Species == presentinall$Genes, ]
					## write name lists to later extract
			for (i in 4:ncol(onetoone)) {
			writeLines(onetoone[,i], paste0(names(onetoone[i]), ".names.txt"))
			}
			dir.create("families_tree")		
			## and gene lists grouped by families
			for (i in 1:nrow(onetoone)) {
			oo <- as.matrix(onetoone)
			writeLines(oo[i, 4:ncol(onetoone)], paste0("families_tree/family", as.character(i), ".names.txt"))
			}
			```
#420 single-copy orhologs present in 36 species. Not that bad...


## now extract the sequences from all fasta files (by species so far)
cd ../
mkdir nt_alignment
for file in *cds; do xargs faidx $file < cdhit_aa/$(basename $file ".cds").pep.repr.faa.names.txt | fasta_formatter >nt_alignment/$file.oto.fa ; done

cd nt_alignment/
cat *oto.fa >all.fasta
## and now extract for each protein group one-by-one
cp ../cdhit_aa/families_tree/family*names.txt .
for names in family*names.txt; do xargs faidx all.fasta < $names > $names.fasta ; done

## align each protein group
for multifasta in family*fasta; do mafft --auto $multifasta >$multifasta.aln; done

##trim alignment and make it single-line
for alignment in *fasta.aln; do trimal -automated1 -in $alignment -out $alignment.trim.aln; done
for tralignment in *trim.aln; do seqkit seq -w 0 $tralignment > $tralignment.sl.aln; done

##concatenate...
##rename each sequence
for file in *sl.aln; do seqkit replace -p ".*" -r '{nr}' --line-width 0 $file >$file.nums.aln; done
seqkit concat *sl.aln.nums.aln >all_concat.aln

#faSize -detailed all_concat.aln #455k, and all have the same length

nano species.tbl

## here how this file should look
cat species.tbl

# 1       Acanthogammarus_godlewskii
# 2       Asprogammarus_rhodophthalmus
# 3       Boeckaxelia_carpenterii
# 4       Brandtia_latissima
# 5       Caprella_sp
# 6       Carinurus_bicarinatus
# 7       Dorogostaiskia_parasitica
# 8       Echinogammarus_berilloni
# 9       Echiuropus_macronychus
# 10      Eulimnogammarus_cyaneus
# 11      Gammarus_chevreuxi
# 12      Gammarus_fossarum
# 13      Gammarus_lacustris
# 14      Gammarus_minus
# 15      Gammarus_pisinnus
# 16      Gammarus_pulex
# 17      Gammarus_wautieri
# 18      Gmelinoides_fasciatus
# 19      Heterogammarus_sophianosii
# 20      Hirondellea_gigas
# 21      Homalogammarus_brandtii
# 22      Hyalellopsis_costata
# 23      Linevichella_vortex
# 24      Macropereiopus_parvus
# 25      Marinogammarus_marinus
# 26      Micruropus_wahlii
# 27      Odontogammarus_calcaratus
# 28      Ommatogammarus_albinus
# 29      Oxyacanthus_curtus
# 30      Pachyschesis_branchialis
# 31      Palicarinus_puzyllii
# 32      Pallasea_cancelloides
# 33      Parapallasea_wosnessenskii
# 34      Parhyale_hawaiensis
# 35      Poekilogammarus_pictoides
# 36      Talitrus_saltator

seqkit replace -p "(.*)" -r '{kv}' --kv-file species.tbl all_concat.aln >all_concat_species.aln

## finally build the tree! 
iqtree -s all_concat_species.aln -abayes -pre all_concat_species_sorted -nt 6 -alrt 1000 -o Caprella_sp


chi-sq failed

#also make a protein-based tree
cd ../
mkdir aa_alignment
for file in *pep; do xargs faidx $file < cdhit_aa/$(basename $file ".pep").pep.repr.faa.names.txt | fasta_formatter >aa_alignment/$file.oto.fa ; done

## and now extract for each protein group one-by-one
cd aa_alignment
cat *oto.fa>all.fasta
cp ../cdhit_aa/families_tree/family*names.txt .
for names in family*names.txt; do xargs faidx all.fasta < $names > $names.fasta ; done

## align each protein group
for multifasta in family*fasta; do mafft --auto $multifasta >$multifasta.aln; done

##trim alignment and make it single-line
for alignment in *fasta.aln; do trimal -automated1 -in $alignment -out $alignment.trim.aln; done
for tralignment in *trim.aln; do seqkit seq -w 0 $tralignment > $tralignment.sl.aln; done

##concatenate...
##rename each sequence
for file in *sl.aln; do seqkit replace -p ".*" -r '{nr}' --line-width 0 $file >$file.nums.aln; done
seqkit concat *sl.aln.nums.aln >all_concat.aln

#faSize -detailed all_concat.aln #455k

nano species.tbl

seqkit replace -p "(.*)" -r '{kv}' --kv-file species.tbl all_concat.aln >all_concat_species.aln

## finally build the tree! 
iqtree -s all_concat_species.aln -abayes -pre all_concat_species_sorted -nt 6 -alrt 1000 -o Caprella_sp

