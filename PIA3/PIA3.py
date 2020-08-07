#!/usr/bin/env python3


import argparse
import argcomplete
import os
import time
import PIA.modules as modules
import ete3


def make_arguments_parser():
    parser = argparse.ArgumentParser(
        prog='Phylogenetically Informed Annotation version 3',
        description='This program is used for phylogenetic analysis')
        
    parser.add_argument('-i', '--input', dest = 'input_file', 
                        help='Input FASTA file OR folder with FASTA files', 
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
    parser.add_argument('-t','--threads', dest = 'threads',  # ADD THIS!!!
                        help='Number of threads',
                        required=False, default=8,
                        type=int)
    parser.add_argument('-aligner', '--aligner_type',
                        default = 'blast', required = False, type=str,
                        help='Use blast or diamond as aligner')
    argcomplete.autocomplete(parser) # for completing by TAB
                                               
    return parser.parse_args()


if __name__ == "__main__":
    
    start = time.time() # start timer
    
    args = make_arguments_parser() # making argparse arguments
    
    paths_list = [] # creating list for programs' paths

    paths_dict = {'transdecoder_orf' : 0, 'transdecoder_cds': 0, 'iqtree' : 0, 'diamond' : 0, 
              'blast' : 0, 'mafft' : 0, 'cdhit' : 0} # crating dictionary for each program (key - program name)

    with open('%s/config/config.ini' % os.path.dirname(os.path.abspath(__file__))) as f: # open config.ini file with program paths
        for line in f:
            if not line.startswith('#') and line != "\n": # if line is not a commented or empty
                paths_list.append(line.rstrip('\n')) # add path to a paths_list

    paths_dict['transdecoder_orf'] = paths_list[0] # path to TransDecoder.LongOrfs script
    paths_dict['transdecoder_cds'] = paths_list[1] # path to TransDecoder.Predict script
    paths_dict['iqtree'] = paths_list[2] # path to IQ-Tree script
    paths_dict['diamond'] = paths_list[3] # path to DIAMOND script
    paths_dict['makeblastdb'] = paths_list[4] # path to makeblastdb script
    paths_dict['blast'] = paths_list[5] # path to blastx script
    paths_dict['mafft'] = paths_list[6] # path to mafft script
    paths_dict['cdhit'] = paths_list[7] # path to cdhit script
    
    # ADD PATH CHECK
    
    # check of argparse arguments
    
    if args.input_file:
        print('Analysis of input FASTA file') # need to check if file exists
    
    if args.all_transcripts:
        transcripts = "all" # assign variable "transcripts" as "all" for function rename_hits
        print('Find all putative transcripts')
    elif args.cds_only:
        transcripts = "cds" # assign variable "transcripts" as "cds" for function rename_hits
        print('Find only putative coding sequences: sequence must start from methionine and has length more or equal of 1/2 of mean database length')

# ------------------------------------------------single mode-----------------------------------------------------------------------------------
        
    if args.single:
        if args.output_folder:
            try:
                os.mkdir(str(args.output_folder)) # try to create folder
            except FileExistsError:
                print("Directory already exists") # if folder exists
            finally:
                os.chdir(args.output_folder) # change working directory to output folder
                path = os.getcwd() # get working directory path
                print("Output directory path: {}".format(path))
    
        file_name = os.path.basename(os.path.normpath(args.input_file)) # get name of analyzed file
        path_out = os.path.join(args.output_folder, file_name) # output path with analyzed file
    
        if os.path.exists(os.path.join(args.output_folder, file_name)): #check if file is already in output folder
            print('File is already in the output folder')
        else:
            print('Copying file to output folder...')
            os.system("cp {} {}".format(args.input_file, args.output_folder)) # if not, copying file in output folder
    
        
        transdecoder_result = modules.run_transdecoder(path_out, paths_dict['transdecoder_orf'], paths_dict['transdecoder_cds']) # run Transdecoder
        
        if args.aligner_type:
            try:
                if args.aligner_type == 'blast': # if blast has been chosen
                    align_result = modules.blast_search(args.database, transdecoder_result, paths_dict['makeblastdb'], paths_dict['blast']) #run blast
                elif args.aligner_type == 'diamond': # if diamond has been chosen
                    align_result = modules.diamond_search(args.database, transdecoder_result, paths_dict['diamond']) #run diamond
            except:
                print('Unknown aligner! Please specify aligner type (blast or diamond)') # if aligner type has been misspelled
                
        cd_hit_result = modules.cd_hit_clust(paths_dict['cdhit']) # run CD-HIT
        translation_result = modules.translate_hits(cd_hit_result) # translate CDS from clustered hits
        renaming_result = modules.rename_hits(file_name,transcripts, args.database) # rename hits according to the filename and filter if CDS mode has been chosen
        initial_phylogeny_result = modules.build_initial_phylogeny(args.database, paths_dict['mafft'], paths_dict['iqtree']) # building initial phylogeny
        mean_dist = modules.calc_mean_dist(initial_phylogeny_result) # calculating absolute mean deviation of branch length
        phylogeny_result = modules.build_phylogeny(args.database, renaming_result, paths_dict['mafft'], paths_dict['iqtree']) # build phyogeny with target hits
        filter_results = modules.filter_distant_seqs(phylogeny_result, mean_dist, renaming_result) # remove branches with lengths more than 4* absolute mean deviation
        tree = ete3.Tree(phylogeny_result) # make Ete3 tree object
        ancestor = tree.search_nodes(name="RHO_Bos_taurus_AAA30674.1")[0] # search for outgroup (CHANGE: USER CAN DEFINE OUTGROUP)
        tree.set_outgroup(ancestor) # set outgroup
        classify = modules.classify_opsins(tree, os.path.splitext(file_name)[0]) #classify opsin by wave length
        write_opsins = modules.write_types(os.path.splitext(file_name)[0]) #write classified opsins into .fasta file
       # write_file = modules.write_into_file(os.path.splitext(file_name)[0]) #write opsins' number into .csv file
        match_seq = modules.match_amino_nucl(os.path.splitext(file_name)[0]) # match amino acid sequence to its nucleotide using file with cds
        os.system('cp query_class_align.fasta.contree {}_class_tree.contree'.format(os.path.splitext(file_name)[0])) # rename final tree
        print('Analysis for {} is done'.format(renaming_result))
    
        print('{}_opsins_class.fasta is a file with amino acid sequences'.format(renaming_result[:-14]))
    
        print('{}_opsins_nt.fasta is a file with nucleotide sequences'.format(renaming_result[:-14]))
    
        
        
# ------------------------------------------------multiple mode-----------------------------------------------------------------------------------        
        
    elif args.multiple: # analyse files from input folder
        for file_name in os.listdir(args.input_file):
            path_in = os.path.abspath(args.input_file) + "/" + file_name
            if args.output_folder:
                try:
                    os.mkdir(str(args.output_folder))
                except FileExistsError:
                    print("Directory already exists")
                finally:
                    os.chdir(args.output_folder)
                    path = os.getcwd()
                    print("Output directory path: {}".format(path))
    
            path_out = args.output_folder+"/"+file_name
    
            if os.path.exists(str(args.output_folder+"/"+file_name)):
                print('File is already in the output folder')
            else:
                print('Copying file to output folder...')
                os.system("cp {} {}".format(path_in, args.output_folder))
    
        
            transdecoder_result = modules.run_transdecoder(path_out, paths_dict['transdecoder_orf'], paths_dict['transdecoder_cds'])
            blast_result = modules.blast_search(args.database, transdecoder_result, paths_dict['diamond'])
            cd_hit_result = modules.cd_hit_clust(paths_dict['cdhit'])
            translation_result = modules.translate_hits(cd_hit_result)
            renaming_result = modules.rename_hits(file_name,transcripts, args.database)
            initial_phylogeny_result = modules.build_initial_phylogeny(args.database, paths_dict['mafft'], paths_dict['iqtree'])
            median_dist = modules.calc_median_dist(initial_phylogeny_result)
            phylogeny_result = modules.build_phylogeny(args.database, renaming_result, paths_dict['mafft'], paths_dict['iqtree'])
            filter_results = modules.filter_distant_seqs(phylogeny_result, median_dist, renaming_result)
     
            tree = ete3.Tree(phylogeny_result)
            ancestor = tree.search_nodes(name="RHO_Bos_taurus_AAA30674.1")[0]
            tree.set_outgroup(ancestor)
            classify = modules.classify_opsins(tree, file_name[:-6])
            write_opsins = modules.write_types(file_name[:-6])
            write_file = modules.write_into_file(file_name[:-6])
            match_seq = modules.match_amino_nucl(file_name[:-6])
            os.system('cp query_class_align.fasta.contree {}_class_tree.contree'.format(file_name[:-6]))
            print('Analysis for {} is done'.format(renaming_result))
    
            print('{}_opsins_class.fasta is a file with amino acid sequences'.format(renaming_result[:-14]))
    
            print('{}_opsins_nt.fasta is a file with nucleotide sequences'.format(renaming_result[:-14]))
    
    if args.delete_intermediate:
        print('Delete intermediate files') # remove unnecessary files
        os.system('rm query_class*')
        os.system('rm pipeliner.*')
        os.system('rm blast*')
        os.system('rm class_align.fasta*')
        os.system('rm -r *transdecoder*')
        os.system('rm user_database*')
        os.system('rm PIA_results.fasta')
            
    end = time.time() # stop timer
    print('Time for analysis:', end - start) # print result time
        
        