# Opsin diversity in amphipod transcriptomes

Companion scripts to the manuscript ... submitted to ... 

## Content

* PIA3
* `other_scripts`: other useful pieces of code for the analysis.
* `results`: reproducible code for figures.

### Identification of opsin genes of Baikal amphipodes

### Environment settings

Snakemake is required for pipeline running. You can

1) Create and activate new conda environment named `smk`

```commandline
conda env create -n smk -c bioconda snakemake=5.20.1
conda activate smk
```

2) Install snakemake in your activated environment

```commandline
conda install -c bioconda snakemake=5.20.1
```

### Run pipeline

```commandline
snakemake -j 8 --use-conda --conda-prefix /path/to/new/conda --config in_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_data out_dir=/path/diversity_of_opsins_in_amphipods/PIA3/test_out db=/path/diversity_of_opsins_in_amphipods/PIA3/classification_opsins_full_aa.fasta cds=True del=True
```