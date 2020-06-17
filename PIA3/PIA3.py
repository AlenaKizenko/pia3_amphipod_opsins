import argparse
import argcomplete
import os
import time
import PIA.modules as modules
import PIA.nucleotide as nucleotide
import PIA.classification as classification
import ete3


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


if __name__ == "__main__":
    
    start = time.time()
    
    args = make_arguments_parser()
    
    # ADD PATH CHECK
    
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
    
    file_name = os.path.basename(os.path.normpath(args.input_file))
    path_out = args.output_folder+"/"+file_name
    
    if os.path.exists(str(args.output_folder+"/"+file_name)):
        print('File is already in the output folder')
    else:
        print('Copying file to output folder...')
        os.system("cp {} {}".format(args.input_file, args.output_folder))
         
    paths_dict = {'transdecoder' : '', 'iqtree' : '', 'diamond' : '', 
              'blast' : '', 'mafft' : '', 'cdhit' : ''}
    paths_list = []

   
    paths_dict = {'transdecoder' : 0, 'iqtree' : 0, 'diamond' : 0, 
              'blast' : 0, 'mafft' : 0, 'cdhit' : 0}
    paths_list = []

   
    with open('%s/config/config.ini' % os.path.dirname(os.path.realpath(__file__))) as f:
        for line in f:
            if not line.startswith('#') and line != "\n":
                paths_list.append(line.rstrip('\n'))

    cnt = 0
    for key in paths_dict.keys():
        paths_dict[key] = paths_list[cnt]
        cnt += 1
        
    transdecoder_result = modules.run_transdecoder(path_out, paths_dict['transdecoder'])
    blast_result = modules.blast_search(args.database, transdecoder_result, paths_dict['diamond'])
    cd_hit_result = modules.cd_hit_clust(paths_dict['cdhit'])
    translation_result = modules.translate_hits(cd_hit_result)
    renaming_result = modules.rename_hits(file_name, args.database)
    initial_phylogeny_result = modules.build_initial_phylogeny(args.database, paths_dict['mafft'], paths_dict['iqtree'])
    median_dist = modules.calc_median_dist(initial_phylogeny_result)
    phylogeny_result = modules.build_phylogeny(args.database, renaming_result, paths_dict['mafft'], paths_dict['iqtree'])
    filter_results = modules.filter_distant_seqs(phylogeny_result, median_dist, renaming_result)
     
    tree = ete3.Tree(phylogeny_result)
    ancestor = tree.search_nodes(name="RHO_Bos_taurus_AAA30674.1")[0]
    tree.set_outgroup(ancestor)
    classify = classification.classify_opsins(tree, file_name[:-6])
    write_opsins = classification.write_types(file_name[:-6])
    write_file = classification.write_into_file(file_name[:-6])
    match_seq = nucleotide.match_amino_nucl(file_name[:-6])
    
    print('Analysis is done')
    print('{}opsins_class.fasta is a file with amino acid sequences'.format(renaming_result))
    print('{}opsins_nt.fasta is a file with nucleotide sequences'.format(renaming_result))
    end = time.time()
    print('Time for analysis:', end - start)

