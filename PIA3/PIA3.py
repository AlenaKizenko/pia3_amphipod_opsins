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
    argcomplete.autocomplete(parser) 
                                               
    return parser.parse_args()


if __name__ == "__main__":
    
    start = time.time()
    
    args = make_arguments_parser()
    
    paths_list = []

    paths_dict = {'transdecoder_orf' : 0, 'transdecoder_cds': 0, 'iqtree' : 0, 'diamond' : 0, 
              'blast' : 0, 'mafft' : 0, 'cdhit' : 0}
    paths_list = []

    with open('%s/config/config.ini' % os.path.dirname(os.path.abspath(__file__))) as f:
        for line in f:
            if not line.startswith('#') and line != "\n":
                paths_list.append(line.rstrip('\n'))

    paths_dict['transdecoder_orf'] = paths_list[0]
    paths_dict['transdecoder_cds'] = paths_list[1]
    paths_dict['iqtree'] = paths_list[2]
    paths_dict['diamond'] = paths_list[3]
    paths_dict['blast'] = paths_list[4]
    paths_dict['mafft'] = paths_list[5]
    paths_dict['cdhit'] = paths_list[6]
    
    # ADD PATH CHECK
    
    if args.input_file:
        print('Analysis of input FASTA file')
    
    if args.all_transcripts:
        transcripts = "all"
        print('Find all putative transcripts')
    elif args.cds_only:
        transcripts = "cds"
        print('Find only putative coding sequences: sequence must start from methionine and has length more or equal of 1/2 of mean database length')
        
    if args.single:
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
        classify = classification.classify_opsins(tree, file_name[:-6])
        write_opsins = classification.write_types(file_name[:-6])
        write_file = classification.write_into_file(file_name[:-6])
        match_seq = nucleotide.match_amino_nucl(file_name[:-6])
        os.system('cp query_class_align.fasta.contree {}_class_tree.contree'.format(file_name[:-6]))
        print('Analysis for {} is done'.format(renaming_result))
    
        print('{}_opsins_class.fasta is a file with amino acid sequences'.format(renaming_result[:-14]))
    
        print('{}_opsins_nt.fasta is a file with nucleotide sequences'.format(renaming_result[:-14]))
    
        if args.delete_intermediate:
            print('Delete intermediate files')
            os.system('rm query_class*')
            os.system('rm pipeliner.*')
            os.system('rm blast*')
            os.system('rm class_align.fasta*')
            os.system('rm -r *transdecoder*')
            os.system('rm user_database.dmnd')
            os.system('rm PIA_results.fasta')
        end = time.time()
        print('Time for analysis:', end - start)
        
    elif args.multiple:
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
            classify = classification.classify_opsins(tree, file_name[:-6])
            write_opsins = classification.write_types(file_name[:-6])
            write_file = classification.write_into_file(file_name[:-6])
            match_seq = nucleotide.match_amino_nucl(file_name[:-6])
            os.system('cp query_class_align.fasta.contree {}_class_tree.contree'.format(file_name[:-6]))
            print('Analysis for {} is done'.format(renaming_result))
    
            print('{}_opsins_class.fasta is a file with amino acid sequences'.format(renaming_result[:-14]))
    
            print('{}_opsins_nt.fasta is a file with nucleotide sequences'.format(renaming_result[:-14]))
    
            if args.delete_intermediate:
                print('Delete intermediate files')
                os.system('rm query_class*')
                os.system('rm pipeliner.*')
                os.system('rm blast*')
                os.system('rm class_align.fasta*')
                os.system('rm -r *transdecoder*')
                os.system('rm user_database.dmnd')
                os.system('rm PIA_results.fasta')
            
            end = time.time()
            print('Time for analysis:', end - start)
        
        