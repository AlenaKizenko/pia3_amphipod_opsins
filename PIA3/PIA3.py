#!/usr/bin/env python3


import argparse
import os
import time

import PIA.modules as modules
import ete3


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
    group_seqs = parser.add_mutually_exclusive_group(required=True)
    group_seqs.add_argument('-all', '--all_transcripts',
                            help='Perform BLAST search on all transcripts',
                            action="store_true")
    group_seqs.add_argument('-cds', '--cds_only',
                            help='Perform BLAST search only on coding sequences',
                            action="store_true")
    parser.add_argument('-del', '--delete_intermediate',
                        action='store_true', help='Delete intermediate files')
    parser.set_defaults(feature=True)
    parser.add_argument('-o', '--output', dest='output_folder',
                        help='Output folder name',
                        required=True,
                        type=str)
    parser.add_argument('-t', '--threads', dest='threads',
                        help='Number of threads',
                        required=False, default=8,
                        type=int)
    parser.add_argument('-cd_h', '--cd_hit', dest='cd_hit',
                        help='Clustering identity threshold',
                        required=False, default=0.95,
                        type=float)
    parser.add_argument('-aligner', '--aligner_type',
                        default='blast', required=False, type=str,
                        help='Use blast or diamond as aligner')
    parser.add_argument('-in_phylo', '--initial_phylo', required=True, type=str,
                        help='Initial phylogeny based on input database')
    return parser.parse_args()


if __name__ == "__main__":

    start = time.time()  # start timer

    args = make_arguments_parser()  # making argparse arguments
    # check of argparse arguments

    if args.all_transcripts:
        transcripts = "all"  # assign variable "transcripts" as "all" for function rename_hits
        print('Find all putative transcripts')
    elif args.cds_only:
        transcripts = "cds"  # assign variable "transcripts" as "cds" for function rename_hits
        print(
            'Find only putative coding sequences:'
            'sequence must start from methionine and has length more or equal of 1/2 of mean database length')

    try:
        os.mkdir(args.output_folder)  # try to create folder
    except FileExistsError:
        print(f"{args.output_folder} : directory already exists")  # if folder exists
    finally:
        os.chdir(args.output_folder)  # change working directory to output folder
        path = os.getcwd()  # get working directory path
        print("Output directory path: {}".format(path))

    file_name = os.path.basename(os.path.normpath(args.input_file))
    transdecoder_result = modules.run_transdecoder(args.input_file)  # run Transdecoder
    if args.aligner_type == 'blast':  # if blast has been chosen
        align_result = modules.blast_search(args.database, transdecoder_result, args.threads)  # run blast
    elif args.aligner_type == 'diamond':  # if diamond has been chosen
        align_result = modules.diamond_search(args.database, transdecoder_result, args.threads)  # run diamond
    else:
        assert False, "Unknown aligner! " \
                              "Please specify aligner type (blast or diamond)"  # if aligner type has been misspelled
    cd_hit_result = modules.cd_hit_clust(c=args.cd_hit, num_threads=args.threads)  # run CD-HIT
    translation_result = modules.translate_hits(cd_hit_result)  # translate CDS from clustered hits
    renaming_result = modules.rename_hits(args.input_file, transcripts,
                                          args.database)  # rename hits according to the filename and filter if CDS mode has been chosen
    mean_dist = modules.calc_mean_dist(args.initial_phylo)  # calculating absolute mean deviation of branch length
    phylogeny_result = modules.build_phylogeny(args.database, renaming_result,
                                               args.threads)  # build phyogeny with target hits
    filter_results = modules.filter_distant_seqs(args.input_file, phylogeny_result, mean_dist,
                                                 renaming_result)  # remove branches with lengths more than 4* absolute mean deviation
    tree = ete3.Tree(phylogeny_result)  # make Ete3 tree object
    ancestor = tree.search_nodes(name="RHO_Bos_taurus_AAA30674.1")[
        0]  # search for outgroup (CHANGE: USER CAN DEFINE OUTGROUP)
    tree.set_outgroup(ancestor)  # set outgroup
    classify = modules.classify_opsins(tree, os.path.splitext(file_name)[0])  # classify opsin by wave length
    write_opsins = modules.write_types(os.path.splitext(file_name)[0])  # write classified opsins into .fasta file
    # write_file = modules.write_into_file(os.path.splitext(file_name)[0]) #write opsins' number into .csv file
    match_seq = modules.match_amino_nucl(
        os.path.splitext(file_name)[0])  # match amino acid sequence to its nucleotide using file with cds
    os.system(
        f'mv query_class_align.fasta.contree {os.path.splitext(file_name)[0]}_class_tree.contree')  # rename final tree

    print(f'Analysis for {renaming_result} is done')
    print(f'{renaming_result[:-14]}_opsins_class.fasta is a file with amino acid sequences')
    print(f'{renaming_result[:-14]}_opsins_nucl.fasta is a file with nucleotide sequences')

    if args.delete_intermediate:
        print('Delete intermediate files')  # remove unnecessary files
        os.system('rm query_class*')
        os.system('rm pipeliner.*')
        os.system('rm blast*')
        os.system('rm class_align.fasta*')
        os.system('rm -r *transdecoder*')
        os.system('rm user_database*')
        os.system('rm PIA_results.fasta')
    end = time.time()  # stop timer
    print('Time for analysis:', end - start)  # print result time
