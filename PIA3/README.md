## PIA3

Modified from PIA2 (https://github.com/xibalbanus/PIA2).
Only tested in Ubuntu-based Linux systems but should work anywhere else with minimum adjustments in the installation process.

### Installation and environment settings

This pipeline requires `Conda`. If it is not installed on your computer, you need to:

**1.** Download [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (look [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) to decide which is right for you).

**2.** Run bash installation script like 
```commandline
bash Anaconda3-2020.02-Linux-x86_64.sh
``` 
or 
```commandline
bash Miniconda3-latest-Linux-x86_64.sh
```

**3.** Source bashrc to activate conda.

```commandline
source ~/.bashrc
```

**4.**
`Snakemake` is also required for pipeline running. There are two options; use whichever works for you:

**A:** Create and activate new conda environment named `smk`

```commandline
conda create -y -n smk snakemake=5.21.0 -c bioconda -c conda-forge
conda activate smk
```

or 

**B:** Install conda environment named smk from the `smk.yml` file

```commandline
conda env create -n smk --file smk.yml
conda activate smk
```

It takes some time, but 

### Testing

We recommend that you test your PIA3 installation

**1.** Run PIA3 on test data

```commandline
cd PIA3

snakemake -j 8 --use-conda --conda-prefix CONDA_PREFIX --config in_dir=/path/to/PIA3/test_data out_dir=/path/to/PIA3/test_out/transcriptome db=/path/to/PIA3/test_data/classification_opsins_full_aa.fasta cds=True del=False opsin=True model=LG+F+G4 outgroup=RHO_Bos_taurus_AAA30674.1
```
**2.** Run unit test

```commandline
 python -m unittest test_PIA3.py
```

Intended for own use. Please feel free to use, reuse, modify and contact us if you need help.


### Run the pipeline

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

* `cdhit`: CH-HIT clustering treshhold (if choose 1, CH-HIT clusters only identical sequences)**default 0.95**

 

### Understanding output files

* initial tree built on the db file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/class_align.fasta.contree`

* results correspond to `test_data/header.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/header`

* results correspond to `test_data/Parhyale_hawaiensis_test.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/Parhyale_hawaiensis_test`

* conda environemnt with all required packages: `/path/to/new/conda`


