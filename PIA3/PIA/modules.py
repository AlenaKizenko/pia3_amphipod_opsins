import os
import re
import statistics

import ete3
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord


# import csv


def run_transdecoder(file):
    os.system(f"TransDecoder.LongOrfs -t {file}")  # run TransDecoder.LongOrfs script
    os.system(f"TransDecoder.Predict -t {file}")  # run TransDecoder.Predict script
    base_name = os.path.basename(os.path.normpath(file))
    transdecoder_cds = f'{base_name}.transdecoder.cds'  # define the name of Transdecoder's result file with cds
    return transdecoder_cds  # return name of Transdecoder's result file with cds


def blast_search(db, cds, num_threads=16, e_val=1e-20):
    os.system(f'makeblastdb -in {db} -out user_database -parse_seqids -dbtype prot')  # make blast database
    os.system(
        f'blastx -query {cds} -db user_database -num_threads {num_threads} -outfmt 6 '
        f'-out blast_file.tmp -evalue {e_val}')  # perform blast search
    hits = set()  # create set for blast hits: since names can be repeated, there is no need to use list
    my_records = []  # create list for hits to write them to .fasta file
    with open(f'blast_file.tmp') as blast_hits:  # open blast file
        for blast_hit in blast_hits:
            hits.add(blast_hit.split()[
                         0])  # write blast hits' name to the file (split because of Transdecoder's prolonged names)
    with open(cds) as transcriptome:  # open transcriptome file
        for seq_record in SeqIO.parse(transcriptome, "fasta"):
            for hit in hits:  # iterate on all transcriptes and blast hits
                if str(hit) == str(seq_record.id):  # if names are equal
                    rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # make SeqRecord object
                    my_records.append(rec)  # write seq record to list my_records
        SeqIO.write(my_records, f'blast_hits_nt.fasta', 'fasta')  # write records to .fasta file


def diamond_search(db, cds, num_threads=16, e_val=1e-20):  # the same but for Diamond
    os.system(f'diamond makedb --in {db} -d user_database')
    os.system(f'diamond blastx -q {cds} -d user_database -p {num_threads} -f 6 -o blast_file.tmp -e {e_val}')
    hits = set()
    my_records = []
    with open('blast_file.tmp') as blast_hits:
        for blast_hit in blast_hits:
            hits.add(blast_hit.split()[0])
    with open(cds) as transcriptome:
        for seq_record in SeqIO.parse(transcriptome, "fasta"):
            for hit in hits:
                if str(hit) == str(seq_record.id):
                    rec = SeqRecord(seq_record.seq, seq_record.id, description='')
                    my_records.append(rec)
        SeqIO.write(my_records, 'blast_hits_nt.fasta', 'fasta')


def cd_hit_clust(c=0.95, n=10, d=0, m=1600, num_threads=8):
    print('Performing clustering of found transcripts using CD-hit')
    hits_clust = 'blast_hits_nt_clust.fasta'  # define CD-HIT result file name
    os.system(
        f'cd-hit-est -i blast_hits_nt.fasta -o {hits_clust} '
        f'-c {c} -n {n} -d {d} -M {m} -T {num_threads}')  # run CD-HIT
    return hits_clust  # return D-HIT redult file name


def translate_hits(hits_clust):
    print('Translating clustered hits')
    my_records = []  # create file for translated hits
    hits_res = f'blast_hits_clust_aa.fasta'  # define name of file with translated hits
    with open(hits_clust) as hits_file:  # open file with clustered hits
        for seq_record in SeqIO.parse(hits_file, "fasta"):  # parse file with hits usinf SeqIO
            aa_rec = SeqRecord(seq_record.seq.translate(to_stop=True), seq_record.id,
                               description='')  # translate clustered hit
            my_records.append(aa_rec)  # append list with translated hit
        SeqIO.write(my_records, hits_res, 'fasta')  # write translated hits to file
    return hits_res  # return name of file with translated hits


def rename_hits(species_name, transcripts, db):
    print('Renaming translated hits')
    length_seqs = [len(seq_record2.seq) for seq_record2 in
                   SeqIO.parse(db, "fasta")]  # create list with database sequences' lengths
    mean_length = statistics.mean(length_seqs)  # calculate mean of sequences' lengths
    my_records = []  # create list for writing of renamed hits
    base_name = os.path.splitext(os.path.basename(os.path.normpath(species_name)))[0]
    filename = f'{base_name}_hits_aa.fasta'  # define file name of translated hits (merged with with initial filename)
    for seq_record in SeqIO.parse('blast_hits_clust_aa.fasta', "fasta"):  # parse file with blast hits
        if transcripts == "cds":  # if cds mode has been chosen - rename hits and filter them
            if str(seq_record.seq)[0] == 'M' and len(str(seq_record.seq)) >= \
                    mean_length // 2:  # if seq starts with Met and its len is longer than 1/2 of mean len
                name = seq_record.id[:seq_record.id.find(" ")]  # take part of record name before space
                final = f'{base_name}_{name}'  # define final name of sequence
                rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with sequence records
        elif transcripts == "all":  # if all mode has been chosen - just rename hits
            name = seq_record.id[:seq_record.id.find(" ")]
            print(name)
            final = f'{base_name}_{name}'  # define final name of sequence
            rec = SeqRecord(seq_record.seq, id=final, description='')
            my_records.append(rec)
    SeqIO.write(my_records, filename, 'fasta')
    return filename


def calc_mean_dist(tree):
    print('Calculating median absolute deviation of evolutionary distances')
    tree = ete3.Tree(tree)  # assign Ete3 tree object
    lst = []  # create list for distances
    for leaf in tree.iter_leaves():  # iterate on leaves
        dist = leaf.dist  # assign variable for distance
        lst.append(dist)  # append list with distances
    me = statistics.mean(lst)  # calculate distances mean
    lst_me = []  # create list for abslute deviations from mean
    for i in lst:
        a = abs(i - me)  # from each distance subtract mean and take absolute value
        lst_me.append(a)  # append list with abslute deviation from mean
    return statistics.mean(lst_me) * 4  # return mean absolute deviation * 4


def build_phylogeny(db, filename, num_threads=8):
    print('Building phylogeny')
    os.system(f'cat {db} {filename} > query_class.fasta')  # merge database and hits in one .fatsa file
    os.system(
        f'mafft --thread {num_threads} --inputorder --auto query_class.fasta > '
        f'query_class_align.fasta')  # align database sequences with mafft
    os.system(f'iqtree -s query_class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST')  # build phylogeny with IQ-Tree
    tree = 'query_class_align.fasta.contree'  # define tree name
    return tree  # return tree name


def filter_distant_seqs(input_file, tree_query, dist_dev, query_file):
    print('Filtering distinct hits')
    tree_query = ete3.Tree(tree_query)  # assign Ete3 tree object
    lst_seqs = []  # create list for selected leaves (distance less than 4*absolute mean deviation)
    my_records = []  # create list for
    file_prefix = os.path.splitext(os.path.basename(os.path.normpath(input_file)))[0]  # define file prefix (without _hits_aa.fasta)
    for leaf in tree_query.iter_leaves():  # iterate on leaves of query tree
        if file_prefix in str(
                leaf.name) and leaf.dist < dist_dev:  # if seq name from file == seq name from tree and distance is OK
            lst_seqs.append(str(leaf.name))  # append lst_seqs by name of selected leaf
    with open(query_file) as hits:
        for seq_record in SeqIO.parse(hits, "fasta"):  # parse .fasta file with hits
            for seq_name in lst_seqs:
                seq_record.id = seq_record.id.replace(":", "_")  # replacement because of IQ-Tree specificity
                seq_record.id = seq_record.id.replace("|", "_")  # replacement because of IQ-Tree specificity
                if seq_record.id in seq_name:  # if seq name from file == selected name
                    rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # create SeqRecord object
                    my_records.append(rec)  # append list with seq records
        SeqIO.write(my_records, 'PIA_results.fasta', 'fasta')  # write seq records to file


def check_lysine(alignment, n=296, ref_seq_name='RHO_Bos_taurus_AAA30674.1'):  # Polina's function
    # calculate position in the alignment
    with open(alignment) as aln:
        for seq_record in SeqIO.parse(aln, "fasta"):
            if str(seq_record.id) == ref_seq_name:  # if the ref sequence found
                sequence_str = seq_record.seq
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

    not_opsins = []
    for seq_record in SeqIO.parse(aln, "fasta"):
        print(seq_record.id)
        print(str(seq_record.seq)[number_position])
        if str(seq_record.seq)[number_position] != 'K':
            not_opsins.append(seq_record.id)

    print(not_opsins)


def match_amino_nucl(filename):
    print('Matching protein sequence to nucleotide CDS')
    file_opsins_set = {title for title, seq in SimpleFastaParser(
        open(filename + '_opsins_class.fasta'))}  # create set from seq names from classified opsins' file
    my_records = []  # create list for
    cnt = 1  # count
    with open('blast_hits_nt_clust.fasta') as transcriptome:  # open file with clustered blast hits
        for seq_record in SeqIO.parse(transcriptome, "fasta"):  # parse file with clustered blast hits
            for opsin in file_opsins_set:  # iterate on names from classified opsins' file
                seq_record_id = str(seq_record.id).replace("|", "_")  # replacement because of IQ-Tree specificity
                seq_record_id = seq_record_id.replace(":", "_")  # replacement because of IQ-Tree specificity
                if seq_record_id.split()[0][:-1] in opsin:  # if seq names are equal
                    name = opsin + str(cnt)  # cnt to avoid similar names
                    rec = SeqRecord(seq_record.seq, id=name, description='')  # create SeqRecord object
                    my_records.append(rec)  # append list with seq records
                    cnt += 1
        SeqIO.write(my_records, filename + '_opsins_nucl.fasta', 'fasta')  # write seq records to file


def classify_opsins(tree, species):
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
            leaf.name = str(type_opsin) + '_' + str(leaf.name)  # rename leaf name
    return tree


def write_types(filename):
    with open('PIA_results.fasta') as file:  # open file with final hits
        my_records = []  # create list for
        path_out = str(filename + '_opsins_class.fasta')  # create filename for .fasta file with classified opsins
        for seq_record in SeqIO.parse(file, "fasta"):  # parse file with final hits
            for key, value in opsins_class.items():  # open dict with info about types
                if key in seq_record.id:  # if sequence name in dict == sequence name in file with final hits
                    name = seq_record.id  # define new sequence name
                    final = f'{value}_{name}'  # concatenate opsin type and sequence name
                    rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                    my_records.append(rec)  # append list with seq records
        SeqIO.write(my_records, path_out, 'fasta')  # write classified sequences to file
