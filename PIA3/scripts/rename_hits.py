from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import statistics
import sys
import os


def rename_hits(db, input_file, output_file, species, transcripts):
    seqs_length = [len(seq_record.seq) for seq_record in
                   SeqIO.parse(db, "fasta")]  # create list with sequence lengths of database sequences
    mean_length = statistics.mean(seqs_length)  # calculate mean of sequences lengths
    my_records = []  # create list for writing of renamed hits
    for seq_record in SeqIO.parse(input_file, "fasta"):  # parse file with blast hits
        if transcripts == "cds":  # if cds mode has been chosen - rename hits and filter them
            if str(seq_record.seq)[0] == 'M' and len(str(seq_record.seq)) >= \
                    mean_length // 2:  # if seq starts with Met and its len is longer than 1/2 of mean len
                name = seq_record.id[:seq_record.id.find(" ")]  # take part of record name before space
                final = f'{species}_{name}'  # define final name of sequence
                rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with sequence records
        elif transcripts == "all":  # if all mode has been chosen - just rename hits
            name = seq_record.id[:seq_record.id.find(" ")]
            final = f'{species}_{name}'  # define final name of sequence
            rec = SeqRecord(seq_record.seq, id=final, description='')
            my_records.append(rec)
    SeqIO.write(my_records, output_file, 'fasta')
    #print(my_records)

db = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]
species = sys.argv[4]
mode = sys.argv[5]

renaming_result = rename_hits(db, input_file, output_file, species, mode)

