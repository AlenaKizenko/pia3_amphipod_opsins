import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import ete3
import statistics

# from Bio.SeqIO.FastaIO import SimpleFastaParser


## hard coded paths. I don't care here.
paths_dict = dict()
paths_dict['mafft'] = 'mafft'
paths_dict['iqtree'] = '/home/drozdovapb/lib/iqtree-1.6.12-Linux/bin/iqtree'


## build_phylogeny from pia3

def align(db, filename, path1):
    print('Aligning sequences')
    os.system('cat {} {} > query_class.fasta'.format(db, filename))
    os.system('{} --thread 8 --inputorder --auto query_class.fasta > query_class_align.fasta'.format(path1))
    aln = 'query_class_align.fasta'
    return aln


## the second part of build_phylogeny
def build_phylogeny(db, filename, path2):
    print('Building phylogeny')
    os.system('cat {} {} > query_class.fasta'.format(db, filename))
    os.system('{} -s query_class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST'.format(path2))
    tree = 'query_class_align.fasta.contree'
    return tree


## 
# phylogeny_result = modules.build_phylogeny(args.database, renaming_result, paths_dict['mafft'], paths_dict['iqtree'])

alignment = align('~/lib/PIA3_v20200628/classification_opsins_full_aa.fasta',
                  '~/Research/Projects/Pigments/plots_tables/2020-07-17-extract-from-mult-alignment/G_sp_Ecy_Gfa_Haz_Pha.fa',
                  paths_dict['mafft'])
# phylogeny_result = build_phylogeny('~/lib/PIA3_v20200628/classification_opsins_full_aa.fasta', '~/Research/Projects/Pigments/plots_tables/2020-07-17-extract-from-mult-alignment/G_sp_Ecy_Gfa_Haz_Pha.fa', paths_dict['iqtree'])

print(alignment)


def get_amino_acid(n, ref_seq_name):
    ## calculate position
    with open(alignment) as aln:
        for seq_record in SeqIO.parse(aln, "fasta"):
            if str(seq_record.id) == ref_seq_name:  ## if the ref sequence found
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

    with open(alignment) as aln, open('output' + str(n) +'.txt', 'w') as outtable:
        for seq_record in SeqIO.parse(aln, "fasta"):
            print(seq_record.id)
            ## outtable.write(str(seq_record.id) + '\n')
            print(str(seq_record.seq)[number_position])
            outtable.write(str(seq_record.seq[number_position]) + '\n')


n = 296
ref_seq_name = 'RHO_Bos_taurus_AAA30674.1'

get_amino_acid(n, ref_seq_name)

get_amino_acid(90, ref_seq_name)
get_amino_acid(113, ref_seq_name)
get_amino_acid(118, ref_seq_name)
get_amino_acid(122, ref_seq_name)
get_amino_acid(126, ref_seq_name)
get_amino_acid(164, ref_seq_name)
get_amino_acid(211, ref_seq_name)
get_amino_acid(269, ref_seq_name)
get_amino_acid(292, ref_seq_name)
get_amino_acid(295, ref_seq_name)
get_amino_acid(296, ref_seq_name)
