#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import time
import PIA.modules as modules
import ete3
from Bio import SeqIO
import sys


def make_arguments_parser():
    parser = argparse.ArgumentParser(
        prog='Phylogenetically Informed Annotation version 3',
        description='This program is used for phylogenetic analysis')

    parser.add_argument('-i', '--input', dest='input_file',
                        help='Input FASTA file OR folder with FASTA files',
                        required=True,
                        type=str)
    parser.add_argument('-db', '--database', dest='database',
                        help='Path to database for blast search and EPA',
                        required=True,
                        type=str)
    parser.add_argument('-cds', '--cds_only',
                            help='Perform BLAST search only on coding sequences',
                            action="store_true")
    parser.add_argument('-del', '--delete_intermediate',
                        action='store_true',
                        help='Delete intermediate files')
    parser.set_defaults(feature=True)
    parser.add_argument('-o', '--output', dest='output_folder',
                        help='Output folder name',
                        required=True,
                        type=str)
    parser.add_argument('-t', '--threads', dest='threads',
                        help='Number of threads',
                        required=False, default=8,
                        type=int)
    parser.add_argument('-aligner', '--aligner_type',
                        required=False, type=str,
                        help='Use BLAST or DIAMOND for database search')
    parser.add_argument('-in_phylo', '--initial_phylo',
                        required=True, type=str,
                        help='Initial phylogeny based on input database')
    parser.add_argument('-model', '--model_iqtree',
                        required=False, type=str,
                        help='Model for IQ-Tree maximum likelihood tree building (if known)')
    parser.add_argument('-outgroup', '--tree_outgroup',
                        required=False, type=str,
                        help='Outgroup for phylogenetic tree building; if not defined by user, first sequence from database FASTA file is taken')
    parser.add_argument('-opsin', '--opsins_search',
                        required=False, action='store_true',
                        help='Searching for opsin sequences (MWS, LWS, UV, Vertebrate-like)')
    parser.add_argument('-cdhit', '--cdhit_threshold',
                        required=False, default=0.95,
                        help='CD-HIT similarity rate threshold for merging transcripts')
    return parser.parse_args()


if __name__ == "__main__":

    start = time.time()  # start timer
    
    args = make_arguments_parser()  # making argparse arguments
    # check of argparse arguments

    if args.cds_only:
        transcripts = "cds"  # assign variable "transcripts" as "cds" for function rename_hits
        print(
            'Find only putative coding sequences:'
            'sequence must start from methionine and has length more or equal of 1/2 of mean database length')
    else:
        transcripts = "all"  # assign variable "transcripts" as "all" for function rename_hits
        print('Find all putative transcripts')

    try:
        os.mkdir(args.output_folder)  # try to create folder
    except FileExistsError:
        print(f"{args.output_folder} : directory already exists")  # if folder exists
    finally:
        os.chdir(args.output_folder)  # change working directory to output folder
        path = os.getcwd()  # get working directory path
        print("Output directory path: {}".format(path))
    
    used_command = " ".join(sys.argv)
    file_w_path = open(f'{args.output_folder}/reproduce_command.txt', "w")
    file_w_path.write(used_command)
    file_w_path.close()

    file_name = os.path.basename(os.path.normpath(args.input_file))
    base_name = os.path.splitext(file_name)[0]
    transdecoder_result = modules.run_transdecoder(args.input_file, args.output_folder)  # run Transdecoder
    if args.aligner_type == 'blast':  # if blast has been chosen
        align_result = modules.blast_search(base_name, args.output_folder, args.database, transdecoder_result, args.threads)  # run blast
    elif args.aligner_type == 'diamond':  # if diamond has been chosen
        align_result = modules.diamond_search(base_name, args.output_folder, args.database, transdecoder_result, args.threads)  # run diamond
    else:
        assert False, "Unknown aligner! " \
                              "Please specify aligner type (blast or diamond)"  # if aligner type has been misspelled
    cd_hit_result = modules.cd_hit_clust(base_name, out_dir=args.output_folder, c = args.cdhit, num_threads=args.threads)  # run CD-HIT
    renaming_result = modules.rename_hits(base_name, args.output_folder,
                                          transcripts, args.database, base_name)  # rename hits according to the filename and filter if CDS mode has been chosen
    mean_dist = modules.calc_mean_dist(args.initial_phylo)  # calculating absolute mean deviation of branch length
    
    alignment_result = modules.mafft_align(args.output_folder, args.database, renaming_result, args.threads) # align database and target sequences with MAFFT
    
    if args.opsins_search:
        opsins = True
        opsins_check = modules.check_lysine(args.output_folder, 'query_class_align.fasta', base_name)
        phylogeny_result = modules.build_phylogeny(args.output_folder, args.database, base_name, args.threads, args.model_iqtree, opsins)  # build phylogeny with target hits        

        filter_results = modules.filter_distant_seqs(base_name, phylogeny_result, mean_dist,
                                                 opsins_check)  # remove branches with lengths more than 4* absolute mean deviation
        
        if args.tree_outgroup != 'file':
            outgroup_name = args.tree_outgroup
        else:
            outgroup_name = [seq_record.id for seq_record in
                   SeqIO.parse(args.database, "fasta")][0]

        tree = ete3.Tree(phylogeny_result)  # make Ete3 tree object
        ancestor = tree.search_nodes(name=outgroup_name)[0]  # search for outgroup
        tree.set_outgroup(ancestor)  # set outgroup
        classify = modules.classify_opsins(tree, base_name)  # classify opsin by wave length
        write_opsins = modules.write_types(base_name)  # write classified opsins into .fasta file
       # write_file = modules.write_into_file(base_name) #write opsins' number into .csv file
    
    else:
        opsins = False
        phylogeny_result = modules.build_phylogeny(args.output_folder,  args.database, renaming_result, args.threads, args.model_iqtree)  # build phylogeny with target hits        
        
        filter_results = modules.filter_distant_seqs(base_name, phylogeny_result, mean_dist,
                                                 renaming_result)
    
    match_seq = modules.match_amino_nucl(file_name, opsins)  # match amino acid sequence to its nucleotide sequence using file with cds
    match_con = modules.match_amino_contig(file_name, args.input_file, opsins) # match amino acid sequence to its contig using file with cds
    os.system(
        f'mv query_class_align.fasta.contree {base_name}_class_tree.contree')  # rename final tree

    print(f'Analysis for {base_name} is done')
    print(f'{base_name}_PIA3_aa.fasta is a file with amino acid sequences')
    print(f'{base_name}_PIA3_nucl.fasta is a file with nucleotide sequences')
    print(f'{base_name}_PIA3_contigs.fasta is a file with contig sequences')
    print(f'{base_name}_class_tree.contree is a maximum likelihood phylogenetic tree with target sequences')
    print(f'To reproduce this analysis run command from file reproduce_command.txt')

    if args.delete_intermediate:
        print('Delete intermediate files')  # remove unnecessary files
        os.system('rm query_class*')
        os.system('rm pipeliner.*')
        os.system('rm blast*')
        os.system('rm class_align.fasta*')
        os.system('rm -r *transdecoder*')
        os.system('rm user_database*')
    end = time.time()  # stop timer
    print('Time for analysis:', round((end - start)/60, 2), 'min')  # print result time
