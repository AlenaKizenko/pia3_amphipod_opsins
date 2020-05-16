from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
import sys

filename = sys.argv[2]
path = sys.argv[1]
opsins_file = path + str('/') + filename + '_opsins_class.fasta'
transcriptome_file = path + str('/') + 'blast_hits_nt_clust.fasta'
#transcriptome_file = path + str('/') + filename + '.fasta.transdecoder.cds'
file_opsins_set = {title for title, seq in SimpleFastaParser(open(opsins_file))}
my_records = []
cnt = 1

with open(transcriptome_file) as transcriptome:
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        for opsin in file_opsins_set:
            seq_record_id = str(seq_record.id).replace("|", "_")
            if seq_record_id.split()[0][:-1] in opsin:
                name = opsin+str(cnt)
                rec = SeqRecord(seq_record.seq, id = name, description = '')
                my_records.append(rec)
                cnt += 1
    path_out = str(path+'/' + filename + '_opsins_nucl.fasta')
    SeqIO.write(my_records, path_out, 'fasta')
