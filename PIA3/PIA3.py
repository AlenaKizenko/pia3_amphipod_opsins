import argparse
import argcomplete
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import ete3
import statistics
import time


start = time.time()


def make_arguments_parser():
    parser = argparse.ArgumentParser(
        prog='Phylogenetically Informed Annotation version 3',
        description='This program is used for phylogenetic analysis')
        
    parser.add_argument('-i', '--input', dest = 'input_file', 
                        help='Input FASTA file', 
                        required=True,
                        type=str)
    parser.add_argument('-db', '--database', dest = 'database',
                        help='Path to database for blast search and EPA',
                        required=True,
                        type=str)
    group_mode = parser.add_mutually_exclusive_group(required=True)
    group_mode.add_argument('-s', '--single', 
                        help='Analyse one FASTA file',
                        action="store_true")
    group_mode.add_argument('-m', '--multiple', 
                        help='Analyse multiple FASTA file',
                        action="store_true")
    group_seqs = parser.add_mutually_exclusive_group(required=True)
    group_seqs.add_argument('-all', '--all_transcripts', 
                        help='Perform BLAST search on all transcripts',
                        action="store_true")
    group_seqs.add_argument('-cds', '--cds_only',
                        help='Perform BLAST search only on coding sequences',
                        action="store_true")                       
    parser.add_argument('-del','--delete_intermediate',
                        action='store_true', help = 'Delete intermediate files')
    parser.set_defaults(feature=True)
                        
    parser.add_argument('-o','--output', dest = 'output_folder', 
                        help='Output folder name', 
                        required=True,
                        type=str)
    parser.add_argument('-t','--threads', dest = 'threads',
                        help='Number of threads',
                        required=False, default=8,
                        type=int)
    argcomplete.autocomplete(parser) 
                                               
    return parser.parse_args()
    

def run_transdecoder(file):   
    os.system("TransDecoder.LongOrfs -t {}".format(file))
    os.system("cp {}.transdecoder_dir/longest_orfs.cds orfs.fasta".format(args.input_file))
    #os.system("sed '/^$/d' *.faa > {}_orfs.fasta".format(args.input_file))
    transdecoder_cds = 'orfs.fasta'
    return transdecoder_cds

def blast_search(db, cds):
    os.system('/media/secondary/apps/diamond/diamond makedb --in {} -d user_database'.format(db))
    os.system('/media/secondary/apps/diamond/diamond blastx -q {} -d user_database -p 16 -f 6 -o blast_file.tmp -e 0.0000000001'.format(cds))
    hits = []
    my_records = []
    with open('blast_file.tmp') as blast_hits:
        for blast_hit in blast_hits:
            hits.append(blast_hit.split()[0])
        hits = set(hits)
    with open(cds) as transcriptome:
        for seq_record in SeqIO.parse(transcriptome, "fasta"):
            for hit in hits:
                if str(hit) == str(seq_record.id):
                    rec = SeqRecord(seq_record.seq, seq_record.id, description = '')
                    my_records.append(rec)
        SeqIO.write(my_records, 'blast_hits_nt.fasta', 'fasta')


def cd_hit_clust():
    print('Performing clustering of found transcripts using CD-hit')
    os.system('/media/secondary/apps/cdhit/cd-hit-est -i blast_hits_nt.fasta -o blast_hits_nt_clust.fasta -c 0.95 -n 10 -d 0 -M 16000 -T 8')
    hits_clust = 'blast_hits_nt_clust.fasta'
    return hits_clust
    
def translate_hits(hits_clust):
    print('Translating clustered hits')
    my_records = []
    with open(hits_clust) as opsins_file:
        for seq_record in SeqIO.parse(opsins_file, "fasta"):
            aa_rec = SeqRecord(seq_record.seq.translate(to_stop = True), seq_record.id, description = '')
            my_records.append(aa_rec)
        SeqIO.write(my_records, path + '/' + 'blast_hits_clust_aa.fasta', 'fasta')
    hits_res = 'blast_hits_clust_aa.fasta'
    return hits_res

def rename_hits():
    print('Renaming translated hits')
    with open('blast_hits_clust_aa.fasta') as file:
        my_records = []
        path_out = str(path+'/' + str(args.input_file[:-6]) + '_hits_aa.fasta')
        filename = str(args.input_file[:-6]) + '_hits_aa.fasta'
        for seq_record in SeqIO.parse(file, "fasta"):
            name = seq_record.id
            final = str(args.input_file[:-6])+str('_')+name
            rec = SeqRecord(seq_record.seq, id = final, description = '')
            my_records.append(rec)
        SeqIO.write(my_records, path_out, 'fasta')
        return filename
    
    
def build_initial_phylogeny():
    print('Building initial tree')
    os.system('mafft --thread 8 --inputorder --auto /media/tertiary/Alena_Kizenko/PIA3/classification_opsins_full_aa.fasta > class_align.fasta')
    os.system('/media/secondary/apps/iqtree-1.6.10-Linux/bin/iqtree -s class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST')
    tree = 'class_align.fasta.contree'
    return tree
    
def calc_median_dist(tree):
    print('Calculating median absolute deviation of evolutionary distances')
    tree = ete3.Tree(tree)
    lst = []
    for leaf in tree.iter_leaves():
        dist = leaf.dist
        lst.append(dist)
    med = statistics.median(lst)
    lst_med = []
    for i in lst:
        a = abs(i-med)
        lst_med.append(a)
    return(statistics.median(lst_med)*4)

def build_phylogeny(filename):
    print('Building phylogeny')
    os.system('cat /media/tertiary/Alena_Kizenko/PIA3/classification_opsins_full_aa.fasta {} > query_class.fasta'.format(filename))
    os.system('mafft --thread 8 --inputorder --auto query_class.fasta > query_class_align.fasta')
    os.system('/media/secondary/apps/iqtree-1.6.10-Linux/bin/iqtree -s query_class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST')
    tree = 'query_class_align.fasta.contree'
    return tree


def filter_distant_seqs(tree_query, dist_dev, query_file):
    print('Filtering distinct hits')
    tree_query = ete3.Tree(tree_query)
    lst_seqs = []
    my_records = []
    file_prefix = query_file[:-14]
    for leaf in tree_query.iter_leaves():
        if file_prefix in str(leaf.name) and leaf.dist < dist_dev:
            lst_seqs.append(str(leaf.name))
    with open(query_file) as hits:
        for seq_record in SeqIO.parse(hits, "fasta"):
            for seq_name in lst_seqs:
                seq_record.id = seq_record.id.replace(":", "_")
                seq_record.id = seq_record.id.replace("|", "_")
                if seq_record.id in seq_name:
                   # print(seq_record.id)
                    rec = SeqRecord(seq_record.seq, seq_record.id, description = '')
                    my_records.append(rec)
        SeqIO.write(my_records, 'PIA_results.fasta', 'fasta')

    
if __name__ == "__main__":
    
    args = make_arguments_parser()
    
    if args.input_file:
        print('Analysis of input FASTA file')
    
    if args.delete_intermediate:
        print('Delete intermediate files')
    
    if args.all_transcripts:
        print('Find all putative transcripts')
    elif args.cds_only:
        print('Find only putative coding sequences')
    
    if args.output_folder:
        try:
            os.mkdir(str(args.output_folder))
        except FileExistsError:
            print("Directory already exists")
        finally:
            os.chdir(args.output_folder)
            path = os.getcwd()
            print("Output directory path: {}".format(path))

transdecoder_result = run_transdecoder(args.input_file)
blast_result = blast_search(args.database, transdecoder_result)
cd_hit_result = cd_hit_clust()
translation_result = translate_hits(cd_hit_result)
renaming_result = rename_hits()
initial_phylogeny_result = build_initial_phylogeny()
median_dist = calc_median_dist(initial_phylogeny_result)
phylogeny_result = build_phylogeny(renaming_result)
filter_results = filter_distant_seqs(phylogeny_result, median_dist, renaming_result)

end = time.time()
print(end - start)
