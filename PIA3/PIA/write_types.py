def write_types(filename):
    print('Writing prefixes to opsin\'s sequences')
    with open(f'{filename}_opsins.fasta') as file:  # open file with final hits
        my_records = []  # create list for seqs, which passed selection by median abs dev
        for seq_record in SeqIO.parse(file, "fasta"):  # parse file with final hits
            for key, value in opsins_class.items():  # open dict with info about types
                if key in seq_record.id:  # if sequence name in dict == sequence name in file with final hits
                    name = seq_record.id  # define new sequence name
                    final = f'{value}_{name}'  # concatenate opsin type and sequence name
                    rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                    my_records.append(rec)  # append list with seq records
        SeqIO.write(my_records, f'{filename}_PIA3_aa.fasta', 'fasta')  # write classified sequences to file
    os.system(f'rm {filename}_opsins.fasta')
    print('Writing prefixes is completed')