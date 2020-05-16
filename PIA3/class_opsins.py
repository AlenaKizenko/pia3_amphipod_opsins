import ete3
import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

tree_file = sys.argv[1]
filename = sys.argv[2]
path = sys.argv[3]

def classify_opsins(tree, species):
    for leaf in tree.iter_leaves():
        if 'RHO' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name = 'OG'
                leaf = leaf.up
        elif 'MEL' in str(leaf) or 'PER' in str(leaf) or 'TMT' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name  = 'VERL'
                leaf = leaf.up
        elif 'MWS' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name = 'MWS'
                leaf = leaf.up        
        elif 'LWS' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name = 'LWS'
                leaf = leaf.up
        elif 'SWS' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name = 'SWS'
                leaf = leaf.up
        elif 'UV' in str(leaf):
            while leaf.up:
                if str(leaf.name)== '':
                    leaf.name = 'UV'
                leaf = leaf.up
            
    for leaf in tree.iter_leaves():
        if species in str(leaf):
            parent = leaf.up
            while str(parent.name) == '':
                parent = parent.up
            type_opsin = str(parent.name)
            opsins_class[leaf.name] = type_opsin
            leaf.name = str(type_opsin) + '_' + str(leaf.name)
 #   print(tree.get_ascii(show_internal=True))
    return tree


def write_types(opsins_class, path, filename):
    with open(path+'/' + 'PIA_results.fasta') as file:
        my_records = []
        path_out = str(path+'/' + filename + '_opsins_class.fasta')
        for seq_record in SeqIO.parse(file, "fasta"):
            for key, value in opsins_class.items():
                if key in seq_record.id:
                    name = seq_record.id
                    final = str(value)+str('_')+name
                    rec = SeqRecord(seq_record.seq, id = final, description = '')
                    my_records.append(rec)
        SeqIO.write(my_records, path_out, 'fasta')

def write_into_file(opsins_class, opsins_number, path, filename):
    for key, value in opsins_class.items():
        for key1, value1 in opsins_number.items():
            if str(value) == str(key1):
                opsins_number[key1] +=1
    with open(path+'/' + str(filename) + '_opsins_class.csv', 'w') as csvfile:
        fieldnames = ['Opsin type', 'Number of found opsins']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key2, value2 in opsins_number.items():
            writer.writerow({'Opsin type': key2, 'Number of found opsins': value2})
            

opsins_class = {}
opsins_number = {'LWS': 0, 'SWS' : 0, 'UV' : 0, 'VERL': 0, 'OG' : 0, 'MWS' : 0}    
tree = ete3.Tree(tree_file)
ancestor = tree.search_nodes(name="RHO_Bos_taurus_AAA30674.1")[0]
tree.set_outgroup(ancestor)
classify = classify_opsins(tree, filename)
write_opsins = write_types(opsins_class, path, filename)
write_types = write_into_file(opsins_class, opsins_number, path, filename)