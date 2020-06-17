## PIA3

Modified from PIA2 (https://github.com/xibalbanus/PIA2).
Only tested on Linux systems.

### Prerequisites
1. `python3`.
2. Standard python packages: `os`, `argparse`, `argcomplete`, `statistics`, `time`. 
3. Specialized python packages: `Biopython` and `ete3`.
Be careful: it's `import Bio`, but in order to download it, you need to (e.g. `pip3 install biopython`)).
4. `TransDecoder`
5. `diamond`
6. `cd-hit`.


### Installation
1. Get and unpack the `pia3_xx.gzip` archive or get the PIA3 folder.
2. Allow execution of the scripts
    `chmod +x *sh`
    `chmod +x *py`
3. Write paths of all listed dependencies in `config.ini` file.

### Running

* adjust the database if necessary. This one was designed for use with Amphipoda sequences.


Intended for own use. Please feel free to use, reuse, modify and contact us if you need help.
