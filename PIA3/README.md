## PIA3

Modified from PIA2 (https://github.com/xibalbanus/PIA2).
Only tested on Linux systems.

### Prerequisites
1. `python3`;
2. Standard python packages: `os`, `argparse`, `argcomplete`, `statistics`, `time`;
3. Specialized python packages: `Biopython` and `ete3`;
Be careful: it's `import Bio`, but you need to install it as `biopython` (e.g. `pip3 install biopython`)).
4. `TransDecoder`;
5. `diamond`;
*Please note that we got the best results with `diamond` \le 0.9.32. The later versions lose some of the results.
6. `cd-hit`.


### Installation
1. Download and unpack the `pia3_xx.gzip` archive or get the PIA3 folder. Navigate to the PIA3 folder.
2. Allow execution of the scripts
    `chmod +x PIA3_run.sh`
    `chmod +x PIA3.py`
3. Write paths of all listed dependencies in the `config/config.ini` file.

### Running

* Adjust the database if necessary. The one provided was designed for use with Amphipoda sequences.
* Run: python3 PIA3.py -i INPUT_FILE -db DATABASE_FILE (-s | -m) (-all | -cds) [-del] -o OUTPUT_FOLDER [-t THREADS]

   ``` 
   This program is used for phylogenetic analysis
   
   optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input FASTA file OR folder with FASTA files
  -db DATABASE, --database DATABASE
                        Path to database for blast search and EPA
  -opsin, --opsins_search
                        Searching for opsin sequences (MWS, LWS, UV, Vertebrate-like)                      
  -cds, --cds_only      Perform BLAST search only on coding sequences
  -del, --delete_intermediate
                        Delete intermediate files
  -model MODEL, --model_iqtree MODEL
                        Model for IQ-Tree maximum likelihood tree building (if known)
  -aligner ALIGNER, --aligner_type ALIGNER
                        Use blast or diamond as aligner
  -outgroup OUTGROUP, --tree_outgroup OUTGROUP
                        Outgroup for phylogenetic tree building; if not defined by user, first sequence from database FASTA file is taken
  -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                        Output folder name
  -t THREADS, --threads THREADS
                        Number of threads
                        ```

Intended for own use. Please feel free to use, reuse, modify and contact us if you need help.
