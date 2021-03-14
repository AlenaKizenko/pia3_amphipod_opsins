from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser


def match_amino_nucl(filename):
    print('Matching protein sequence to nucleotide CDS')
    file_opsins_set = {title for title, seq in SimpleFastaParser(open(filename + '_opsins_class.fasta'))}
    my_records = []
    cnt = 1
    with open('blast_hits_nt_clust.fasta') as transcriptome:
        for seq_record in SeqIO.parse(transcriptome, "fasta"):
            for opsin in file_opsins_set:
                seq_record_id = str(seq_record.id).replace("|", "_")
                if seq_record_id.split()[0][:-1] in opsin:
                    name = opsin+str(cnt)
                    rec = SeqRecord(seq_record.seq, id = name, description = '')
                    my_records.append(rec)
                    cnt += 1
        SeqIO.write(my_records, filename + '_opsins_nucl.fasta', 'fasta')