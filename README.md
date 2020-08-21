# Opsin diversity in amphipod transcriptomes

Companion scripts to the manuscript ... submitted to ... 

## Content

* PIA3
* `other_scripts`: other useful pieces of code for the analysis.
* `results`: reproducible code for figures.

### Identification of opsin genes of Baikal amphipodes

### Environment settings

Conda is required for pipeline running. If it is not installed on your computer, you need:

**1.** Download [Anaconda](https://www.anaconda.com/products/individual)

**2.** Run bash installation script

```commandline
bash Anaconda3-2020.02-Linux-x86_64.sh
```
**3.** Activate shell

```commandline
source ~/.bashrc
```

Snakemake is required for pipeline running. You can

**A:** Create and activate new conda environment named `smk`

```commandline
conda create -y -n smk snakemake=5.21.0 -c bioconda -c conda-forge
conda activate smk
```

or 

**B:** Install conda environment named smk from `smk.yml` file

```commandline
conda env create -n smk --file smk.yml
conda activate smk
```

### Run pipeline

```commandline
cd PIA3
snakemake -j 8 --use-conda --conda-prefix /path/to/new/conda --config in_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_data out_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_out db=/path/diversity_of_opsins_in_amphipods/PIA3/classification_opsins_full_aa.fasta cds=True del=True
```

* `conda-prefix`: where do you want to install env with all required packages **required**

* `in_dir`: path to folder with input reference `.fasta` file(s) **required**

* `out_dir`: output directory path **required**

* `db`: path to database **required**

* `cds`: perform BLAST search only on coding sequences (longer than 1/2 of mean database sequence and starting fron methionine) **default True**

* `del`: delete intermediate files **default True**

* `aligner`: use BLAST or DIAMOND for database search **default blast**

* `model`: model for IQ-Tree maximum likelihood tree building (if known) **default TEST**

* `outgroup`: outgroup for phylogenetic tree building; if not defined by user, first sequence from database FASTA file is taken

* `opsin`: searching for opsin sequences (MWS, LWS, UV, Vertebrate-like) **default True**
 

### Output

* initial tree built on the db file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/class_align.fasta.contree`

* results correspond to `test_data/header.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/header`

* results correspond to `test_data/Parhyale_hawaiensis_test.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/Parhyale_hawaiensis_test`

* conda environemnt with all required packages: `/path/to/new/conda`
