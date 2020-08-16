# Opsin diversity in amphipod transcriptomes

Companion scripts to the manuscript ... submitted to ... 

## Content

* PIA3
* `other_scripts`: other useful pieces of code for the analysis.
* `results`: reproducible code for figures.

### Identification of opsin genes of Baikal amphipodes

### Environment settings

Snakemake is required for pipeline running. You can

**A:** Create and activate new conda environment named `smk`

```commandline
conda create -y -n smk snakemake -c bioconda
conda activate smk
```

**B:** Install snakemake in your activated environment

```commandline
conda install -y -c bioconda snakemake=5.20.1
```

**C:** Install conda environment named smk from `smk.yml` file

```commandline
conda env create -y -n smk --file smk.yml
```

### Run pipeline

```commandline
cd PIA3
snakemake -j 8 --use-conda --conda-prefix /path/to/new/conda --config in_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_data out_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_out db=/path/diversity_of_opsins_in_amphipods/PIA3/classification_opsins_full_aa.fasta cds=True del=True
```

* `conda-prefix`: where do you want to install env with all required packages

* `in_dir`: path to folder with input reference `.fasta` file(s)

* `out_dir`: output directory path

* `db`: path to database

* other arguments: see `python PIA3/PIA3.py --help`

### Output

* initial tree built on the db file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/class_align.fasta.contree`

* results correspond to `test_data/header.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/header`

* results correspond to `test_data/tail.fasta` file: `/path/diversity_of_opsins_in_amphipods/PIA3/test_out/tail`

* conda environemnt with all required packages: `/path/to/new/conda`