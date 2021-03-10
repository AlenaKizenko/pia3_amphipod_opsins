#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import statistics
from Bio import AlignIO
import ete3
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
# import csv


def run_transdecoder(input_file, out_dir):
    print('Running Transdecoder to define CDS in contigs')
    os.system(f"TransDecoder.LongOrfs -t {input_file}")  # run TransDecoder.LongOrfs script
    os.system(f"TransDecoder.Predict -t {input_file} --single_best_only")  # run TransDecoder.Predict script
    file_name = os.path.basename(os.path.normpath(input_file))
    transdecoder_pep = f'{out_dir}/{file_name}.transdecoder.pep'  # define the name of Transdecoder's result file with cds
    print('Transdecoder run is completed')
    return transdecoder_pep  # return name of Transdecoder's result file with cds


def blast_search(basename, out_dir: str, db, cds, num_threads=16, e_val=1e-10):
    print('Performing BLAST search to define CDS in contigs')
    os.system(f'makeblastdb -in {db} -out user_database -parse_seqids -dbtype prot')  # make blast database
    os.system(
        f'blastp -query {cds} -db user_database -num_threads {num_threads} -outfmt 6 '
        f'-out {out_dir}/{basename}_blast_file.tmp -evalue {e_val}')  # perform blast search
    hits = set()  # create set for blast hits: since names can be repeated, there is no need to use list
    my_records = []  # create list for hits to write them to .fasta file
    with open(f'{out_dir}/{basename}_blast_file.tmp') as blast_hits:  # open blast file
        for blast_hit in blast_hits:
            hits.add(blast_hit.split()[
                         0])  # write blast hits' name to the file (split because of Transdecoder's prolonged names)
    for seq_record in SeqIO.parse(cds, "fasta"):
        for hit in hits:  # iterate on all transcriptes and blast hits
            if str(hit) == str(seq_record.id):  # if names are equal
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # make SeqRecord object
                my_records.append(rec)  # write seq record to list my_records
    SeqIO.write(my_records, f'{out_dir}/{basename}_blast_hits.fasta', 'fasta')  # write records to .fasta file
    print('BLAST search is completed')


def diamond_search(basename, out_dir: str, db, cds, num_threads=16, e_val=1e-10):  # the same but for Diamond
    print('Performing DIAMOND search to define CDS in contigs')
    os.system(f'diamond makedb --in {db} -d user_database')
    os.system(f'diamond blastp -q {cds} -d user_database -p {num_threads} -f 6 -o {out_dir}/{basename}_blast_file.tmp -e {e_val}')
    hits = set()
    my_records = []
    with open(f'{out_dir}/{basename}_blast_file.tmp') as blast_hits:
        for blast_hit in blast_hits:
            hits.add(blast_hit.split()[0])
    for seq_record in SeqIO.parse(cds, "fasta"):
        for hit in hits:
            if str(hit) == str(seq_record.id):
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')
                my_records.append(rec)
    SeqIO.write(my_records, f'{out_dir}/{basename}_blast_hits.fasta', 'fasta')
    print('DIAMOND search is completed')


def cd_hit_clust(basename, out_dir: str, c, n=5, m=2000, num_threads=8):
    print('Performing clustering of found transcripts using CD-hit')
    hits_clust = f'{out_dir}/{basename}_blast_hits_clust.fasta'  # define CD-HIT result file name
    os.system(
        f'cd-hit -i {out_dir}/{basename}_blast_hits.fasta -o {hits_clust} '
        f'-c {c} -n {n} -M {m} -T {num_threads}')  # run CD-HIT
    print('CD-HIT clustering is completed')
    return hits_clust  # return D-HIT redult file name


def rename_hits(species_name, out_dir, transcripts, db, basename):
    print('Renaming hits')
    length_seqs = [len(seq_record2.seq) for seq_record2 in
                   SeqIO.parse(db, "fasta")]  # create list with database sequences' lengths
    mean_length = statistics.mean(length_seqs)  # calculate mean of sequences' lengths
    my_records = []  # create list for writing of renamed hits
    filename = f'{out_dir}/{species_name}_hits.fasta'  # define file name of translated hits (merged with with initial filename)
    for seq_record in SeqIO.parse(f'{out_dir}/{basename}_blast_hits_clust.fasta', "fasta"):  # parse file with blast hits
        if transcripts == "cds":  # if cds mode has been chosen - rename hits and filter them
            if str(seq_record.seq)[0] == 'M' and len(str(seq_record.seq)) >= \
                    mean_length // 2:  # if seq starts with Met and its len is longer than 1/2 of mean len
                name = seq_record.id[:seq_record.id.find(" ")]  # take part of record name before space
                final = f'{species_name}_{name}'  # define final name of sequence
                rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with sequence records
        elif transcripts == "all":  # if all mode has been chosen - just rename hits
            name = seq_record.id[:seq_record.id.find(" ")]
            print(name)
            final = f'{species_name}_{name}'  # define final name of sequence
            rec = SeqRecord(seq_record.seq, id=final, description='')
            my_records.append(rec)
    SeqIO.write(my_records, filename, 'fasta')
    print('Renaming is completed')
    return filename


def calc_mean_dist(tree):
    print('Calculating median absolute deviation of evolutionary distances')
    tree = ete3.Tree(tree)  # assign Ete3 tree object
    lst = []  # create list for distances
    for leaf in tree.iter_leaves():  # iterate on leaves
        dist = leaf.dist  # assign variable for distance
        lst.append(dist)  # append list with distances
    me = statistics.mean(lst)  # calculate distances mean
    lst_me = []  # create list for absolute deviations from mean
    for i in lst:
        a = abs(i - me)  # from each distance subtract mean and take absolute value
        lst_me.append(a)  # append list with absolute deviation from mean
    print('Calculations are completed')
    return statistics.mean(lst_me) * 4  # return mean absolute deviation * 4


def mafft_align(out_dir: str, db, filename, num_threads=8):
    print('Aligning sequences with MAFFT')
    os.system(f'cat {db} {filename} > {out_dir}/query_class.fasta')  # merge database and hits in one .fasta file
    os.system(
        f'mafft --thread {num_threads} --inputorder --auto {out_dir}/query_class.fasta > '
        f'{out_dir}/query_class_align.fasta')  # align database sequences with mafft
    print('Aligning is completed')
    

def build_phylogeny(out_dir: str, db, filename, num_threads=8, model = 'TEST', opsins = False):
    print('Building phylogeny')
    if opsins == False:
        os.system(f'iqtree -s {out_dir}/query_class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m {model}')  # build phylogeny with IQ-Tree
        tree = f'{out_dir}/query_class_align.fasta.contree'  # define tree name
    elif opsins == True:
        os.system(f'cat {db} {filename}_opsins.fasta > {out_dir}/query_class_opsins.fasta')  # merge database and hits in one .fasta file
        os.system(
                f'mafft --thread {num_threads} --inputorder --auto {out_dir}/query_class_opsins.fasta > '
                f'{out_dir}/query_class_opsins_align.fasta')  # align database sequences with mafft
        os.system(f'iqtree -s {out_dir}/query_class_opsins_align.fasta -nt AUTO -t RANDOM -bb 1000 -m {model}')  # build phylogeny with IQ-Tree
        tree = f'{out_dir}/query_class_opsins_align.fasta.contree'  # define tree name
    print('Building of phylogenetic tree is completed')
    return tree  # return tree name


def filter_distant_seqs(input_file, tree_query, dist_dev, query_file):
    print('Filtering distinct hits')
    tree_query = ete3.Tree(tree_query)  # assign Ete3 tree object
    lst_seqs = []  # create list for selected leaves (distance less than 4*absolute mean deviation)
    my_records = []  # create list for selected sequences
    for leaf in tree_query.iter_leaves():  # iterate on leaves of query tree
        if input_file in str(
                leaf.name) and leaf.dist < dist_dev:  # if seq name from file == seq name from tree and distance is OK
            lst_seqs.append(str(leaf.name))  # append lst_seqs by name of selected leaf
    for seq_record in SeqIO.parse(query_file, "fasta"):  # parse .fasta file with hits
        for seq_name in lst_seqs:
            seq_record.id = seq_record.id.replace(":", "_")  # replacement because of IQ-Tree specificity
            seq_record.id = seq_record.id.replace("|", "_")  # replacement because of IQ-Tree specificity
            if seq_record.id in seq_name:  # if seq name from file == selected name
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, f'{input_file}_PIA3_aa.fasta', 'fasta')  # write seq records to file
    print('Filtering is completed')



def check_lysine(out_dir, alignment, filename, ref_seq_name = 'RHO_Bos_taurus_AAA30674.1', n=296): # calculate position in the alignment
    print('Checking lysine position in opsins\' sequences')
    with open(alignment) as aln:
        for record in AlignIO.read(aln, "fasta"):
            if str(record.id) == ref_seq_name:  # if the ref sequence found
                sequence_str = record.seq
                letter_counter = 0
                gap_counter = 0
                for letter in sequence_str:
                    if letter_counter < n - 1:
                        if letter == "-":
                            gap_counter += 1
                        else:
                            letter_counter += 1
                break

    number_position = gap_counter + letter_counter

    all_opsins = []
    with open(alignment) as aln:
        for seq_record in AlignIO.read(alignment, "fasta"):
            if str(seq_record.seq)[number_position] == 'K':
                all_opsins.append(seq_record.id)       
    query_opsins = []
    for seq_record in SeqIO.parse(f'{out_dir}/{filename}_hits.fasta', 'fasta'):
        for opsin in all_opsins:
            if opsin in seq_record.id:
                query_opsins.append(seq_record)
    SeqIO.write(query_opsins, f'{filename}_opsins.fasta', 'fasta')
    print('Checking lysine position is completed')
    return f'{filename}_opsins.fasta'


def classify_opsins(tree, species):
    print('Classifying opsins by sensitivity')
    global opsins_class  # assign global variable
    opsins_class = {}  # create dictionary for opsins' types classification
    for leaf in tree.iter_leaves():  # iterate on leaves
        if 'RHO' in str(leaf):  # if RHO (Outgroup - OG) is in reference database sequence name
            while leaf.up:  # move to the root of the tree
                if not str(leaf.name):  # if node name is empty
                    leaf.name = 'OG'  # name node as OG
                leaf = leaf.up  # move to the root of the tree once
        elif 'MEL' in str(leaf) or 'PER' in str(leaf) or 'TMT' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'VERL'
                leaf = leaf.up
        elif 'MWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'MWS'
                leaf = leaf.up
        elif 'LWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'LWS'
                leaf = leaf.up
        elif 'SWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'SWS'
                leaf = leaf.up
        elif 'UV' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'UV'
                leaf = leaf.up

    for leaf in tree.iter_leaves():  # iterate on leaves
        if species in str(leaf):  # if leaf is a hit
            parent = leaf.up  # parent is a node
            while not str(parent.name):  # while node name is empty
                parent = parent.up  # move to the root of the tree
            type_opsin = str(parent.name)  # define type of opsin as parent node name
            opsins_class[leaf.name] = type_opsin  # write type of hit to dict
            leaf.name = f'{type_opsin}_{leaf.name}'  # rename leaf name
    print('Classification is completed')
    return tree


def write_types(filename):
    print('Writing prefixes to opsin\'s sequences')
    with open(f'{filename}_opsins.fasta') as file:  # open file with final hits
        my_records = []  # create list for seqs, which passed selection by median abs dev
        for seq_record in SeqIO.parse(file, "fasta"):  # parse file with final hits
            for key, value in opsins_class.items():  # open dict with info about types
                if key in seq_record.id:  # if sequence name in dict == sequence name in file with final hits
                    name = seq_record.id  # define new sequence name
                    final = f'{value}_{name}'  # concatenate opsin type and sequence name
                    rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                    my_records.append(rec)  # append list with seq records
        SeqIO.write(my_records, f'{filename}_PIA3_aa.fasta', 'fasta')  # write classified sequences to file
    os.system(f'rm {filename}_opsins.fasta')
    print('Writing prefixes is completed')
        
def match_amino_nucl(filename, opsins = False):
    print('Matching protein sequence to nucleotide CDS')
    file_hits_set = {title for title, seq in SimpleFastaParser(
            open(f'{os.path.splitext(filename)[0]}_PIA3_aa.fasta'))}  # create set from seq names from classified opsins' file
    my_records = []  # create list for
    for seq_record in SeqIO.parse(f'{filename}.transdecoder.cds', "fasta"):  # parse file with clustered blast hits
        for opsin in file_hits_set:  # iterate on names from classified opsins' file
            if seq_record.id.split()[0][:-1] in opsin:  # if seq names are equal
                rec = SeqRecord(seq_record.seq, id=opsin, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, f'{os.path.splitext(filename)[0]}_PIA3_nucl.fasta', 'fasta')  # write seq records to file
    print('Matching protein sequences to nucleotide CDS is completed')

def match_amino_contig(filename, transcriptome, opsins = False):
    print('Matching protein sequence to contig in transcriptome')
    file_hits_set = {title for title, seq in SimpleFastaParser(
            open(f'{os.path.splitext(filename)[0]}_PIA3_aa.fasta'))}  # create set from seq names from classified opsins' file
    my_records = []  # create list for
    for seq_record in SeqIO.parse(transcriptome, "fasta"):  # parse transcriptome file
        for opsin in file_hits_set:  # iterate on names from classified opsins' file
            if seq_record.id.split()[0][:-1] in opsin:  # if seq names are equal
                rec = SeqRecord(seq_record.seq, id=opsin, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, f'{os.path.splitext(filename)[0]}_PIA3_contigs.fasta', 'fasta')  # write seq records to file
    print('Matching protein sequences to contigs is completed')
