import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import ete3
import statistics


def run_transdecoder(file, path1, path2):
    os.system("{} -t {}".format(path1, file))
    os.system("{} -t {}".format(path2, file))
    transdecoder_cds = file + '.transdecoder.cds'
    return transdecoder_cds

def blast_search(db, cds, path):
    os.system('{} makedb --in {} -d user_database'.format(path, db))
    os.system('{} blastx -q {} -d user_database -p 16 -f 6 -o blast_file.tmp -e 0.0000000001'.format(path, cds))
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


def cd_hit_clust(path):
    print('Performing clustering of found transcripts using CD-hit')
    os.system('{} -i blast_hits_nt.fasta -o blast_hits_nt_clust.fasta -c 0.95 -n 10 -d 0 -M 16000 -T 8'.format(path))
    hits_clust = 'blast_hits_nt_clust.fasta'
    return hits_clust
    
def translate_hits(hits_clust):
    print('Translating clustered hits')
    my_records = []
    with open(hits_clust) as opsins_file:
        for seq_record in SeqIO.parse(opsins_file, "fasta"):
            aa_rec = SeqRecord(seq_record.seq.translate(to_stop = True), seq_record.id, description = '')
            my_records.append(aa_rec)
        SeqIO.write(my_records, 'blast_hits_clust_aa.fasta', 'fasta')
    hits_res = 'blast_hits_clust_aa.fasta'
    return hits_res

def rename_hits(species_name, transcripts, db):
    print('Renaming translated hits')
    length_seqs = [len(seq_record2.seq) for seq_record2 in SeqIO.parse(db, "fasta")]
    mean_length = statistics.mean(length_seqs)
    my_records = []
    filename = species_name[:-6] + '_hits_aa.fasta'
    for seq_record in SeqIO.parse('blast_hits_clust_aa.fasta', "fasta"):
        if transcripts == "cds":
            if str(seq_record.seq)[0] == 'M' and len(str(seq_record.seq)) >= mean_length//2:
                name = seq_record.id[:seq_record.id.find(" ")]
                print(name)
                final = species_name[:-6] + '_' + name
                rec = SeqRecord(seq_record.seq, id = final, description = '')
                my_records.append(rec)
            else:
                pass
        elif transcripts == "all":
            name = seq_record.id[:seq_record.id.find(" ")]
            print(name)
            final = species_name[:-6] + '_' + name
            rec = SeqRecord(seq_record.seq, id = final, description = '')
            my_records.append(rec)
    SeqIO.write(my_records, filename, 'fasta')
    return filename
    
    
def build_initial_phylogeny(db, path1, path2):
    print('Building initial tree')
    os.system('{} --thread 8 --inputorder --auto {} > class_align.fasta'.format(path1, db))
    os.system('{} -s class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST'.format(path2))
    tree = 'class_align.fasta.contree'
    return tree
    
def calc_median_dist(tree):
    print('Calculating median absolute deviation of evolutionary distances')
    tree = ete3.Tree(tree)
    lst = []
    for leaf in tree.iter_leaves():
        dist = leaf.dist
        lst.append(dist)
    me = statistics.mean(lst)
    lst_me = []
    for i in lst:
        a = abs(i-me)
        lst_me.append(a)
    return(statistics.mean(lst_me)*4)

def build_phylogeny(db, filename, path1, path2):
    print('Building phylogeny')
    os.system('cat {} {} > query_class.fasta'.format(db, filename))
    os.system('{} --thread 8 --inputorder --auto query_class.fasta > query_class_align.fasta'.format(path1))
    os.system('{} -s query_class_align.fasta -nt AUTO -t RANDOM -bb 1000 -m TEST'.format(path2))
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

