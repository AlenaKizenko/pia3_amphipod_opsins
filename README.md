# Opsins\` diversity in transcriptomes of Baikal endemic amphipods

## Project description

Amphipoda is an order of malacostracan crustaceans with no carapace and generally with laterally compressed bodies. Amphipoda\`s species inhabit different areas and depth of seas and freshwaters.
For orientation in water, amphipods use opsins. Opsins, which belong to the subfamily of G protein-coupled transmembrane receptors and form visual pigments together with retinal chromophores, play key roles in animal photoreception. Opsins may be sensitive to different wavelengths of light, depending on the depth at which the crustacean lives and some other factors such as amphipoda species coloring.
According to recent research there were two independent invasions of amphipods into Baikal. It would be interesting to discover which opsins Baikal amphippda have and how natural selection influences opsins\` types representation. 

## Goals and objectives

The goal of this project is to assess the diversity of opsins in transcriptomes of Baikal amphipods.

1. Download and assess quality of published transcriptome assemblies of 66 species.
2. Apply the published PIA2 pipeline (phylogenetically-informed annotation) to identify the opsin genes;
3. Build trees according to the most interesting sequences;
4. Evaluate the effect of selection in specific sequences and predict the receptor properties of the found opsins.

## Methods

All work was performed on server.

### Identification of opsin genes of Baikal amphipodes
Transcriptomes\` assemblies were downloaded and renamed manually. Quality control was done using BUSCO (Benchmarking Universal Single-Copy Orthologs, default parameters). Result statistics were summarized using `extract_statistics_busco.py` script, then .txt file was converted to .csv file; family names were added manually as a new column in this .csv file and histogramm, depicting taxonomy/quality relations, were plotted using `family_and_quality.R` script. Two rounds of PIA (phylogenetically-informed annotation) were performed: first time with E-value 10e-20, second time with E-value 10e-10.
Results of PIA pipeline - identified opsins - were counted and quality/quantity histograms were plotted using script.
Filtered opsins were used for gblocks processing and further phylogenetic Bayesian tree building using IQTree (parameters: -st AA -m TEST -bb 1100 -abayes -nt AUTO). Final Baeysian tree was visualized and colored using iTOL (Interactive tree of life).

### Article results repeating
Three pairs of reads were downloaded from SRA using `sra_download` bash script. Then reads were processed with Trimmomatic using `trimmomatic.sh` bash script. PIA pipeline was applied with default parameters.

## Results

First, transcriptomes were downloaded and their quality was analysed. We found, that quality of assemblies is rather bad.

We tested if quality of trancriptomes depends on taxonomic position. We plotted histogramm of Missing BUSCOs and revealed, that there is no relation between these parameters.

![alt text](https://github.com/AlenaKizenko/diversity_of_opsins_in_amphipods/blob/master/results/family_stats.jpg)

Then we performed PIA pipeline on Baikal amphipods assemblies with different E-value. After second round of PIA, we filtered identified opsins, by removing species with more than 20% of Missing BUSCOs.

Finally, we built phylogenetic tree using Bayesian method. We identified long-wave sensitive opsins and opsin-like proteins. Probaly, Baikal endemic amphipods lost short-wave and ultra-violet sensitive opsins due to the natural selection.

![alt text](https://github.com/AlenaKizenko/diversity_of_opsins_in_amphipods/blob/master/results/bayes_tree_final.jpg)

Moreover, we almost repeted article results. We identified 4 long-wave opsins, 2 short-wave opsins and 0 opsin-like proteins, whether authors of this pipeline identified 3 long-wave opsins, 2 short-wave opsins and 1 opsin-like protein.


## References
1. Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A Flexible Trimmer for Illumina Sequence Data.” Bioinformatics (Oxford, England) 30(15):2114–20.
2. Naumenko, Sergey A., Maria D. Logacheva, Nina V. Popova, Anna V. Klepikova, Aleksey A. Penin, Georgii A. Bazykin, Anna E. Etingova, Nikolai S. Mugue, Alexey S. Kondrashov, and Lev Y. Yampolsky. 2017. “Transcriptome-Based Phylogeny of Endemic Lake Baikal Amphipod Species Flock: Fast Speciation Accompanied by Frequent Episodes of Positive Selection.” Molecular Ecology 26(2):536–53.
3. Pérez-Moreno, Jorge L., Danielle M. DeLeo, Ferran Palero, and Heather D. Bracken-Grissom. 2018. “Phylogenetic Annotation and Genomic Architecture of Opsin Genes in Crustacea.” Hydrobiologia 825(1):159–75.
4. Simão, Felipe A., Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov. 2015. “BUSCO: Assessing Genome Assembly and Annotation Completeness with Single-Copy Orthologs.” Bioinformatics 31(19):3210–12.
5. Letunic, I. and P. Bork. 2007. “Interactive Tree Of Life (ITOL): An Online Tool for Phylogenetic Tree Display and Annotation.” Bioinformatics 23(1):127–28.

