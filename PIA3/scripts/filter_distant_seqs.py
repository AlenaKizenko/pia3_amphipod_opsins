import sys
from ete3 import Tree
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def check_lysine(alignment, ref_seq_name = 'RHO_Bos_taurus_AAA30674.1', n=296): # calculate position in the alignment
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
        for seq_record in AlignIO.read(aln, "fasta"):
            if str(seq_record.seq)[number_position] == 'K':
                seq_record.id = seq_record.id.replace(":", "_")  # replacement because of IQ-Tree specificity
                seq_record.id = seq_record.id.replace("|", "_")
                all_opsins.append(seq_record.id)       
    return all_opsins


def filter_distant_seqs(tree_query, query_file, dist_dev, input_file, output, opsins):
    tree_query = Tree(tree_query)  # assign Ete3 Tree object
    lst_seqs = []  # create list for selected leaves (distance less than 4*absolute mean deviation)
    my_records = []  # create list for selected sequences
    with open(dist_dev) as dist:
        dist = float(dist.readline())
    for leaf in tree_query.iter_leaves():  # iterate on leaves of query tree
        if input_file in str(
                leaf.name) and leaf.dist < dist:  # if seq name from file == seq name from tree and distance is OK
            lst_seqs.append(str(leaf.name))  # append lst_seqs by name of selected leaf
    if opsins == 'True':
        opsins_lst = check_lysine(alignment, ref_seq_name='RHO_Bos_taurus_AAA30674.1', n=296)
        lst_seqs = set(lst_seqs).intersection(opsins_lst)
    
    for seq_record in SeqIO.parse(query_file, "fasta"):  # parse .fasta file with hits
        for seq_name in lst_seqs:
            seq_record.id = seq_record.id.replace(":", "_")  # replacement because of IQ-Tree specificity
            seq_record.id = seq_record.id.replace("|", "_")  # replacement because of IQ-Tree specificity
            if seq_record.id in seq_name:  # if seq name from file == selected name
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, output, 'fasta')  # write seq records to file

tree = sys.argv[1]
query_file = sys.argv[2]
dist_dev = sys.argv[3]
input_file = sys.argv[4]
output = sys.argv[5]
alignment = sys.argv[6]
opsins = sys.argv[7]

result_filtering = filter_distant_seqs(tree, query_file, dist_dev, input_file, output, opsins)


