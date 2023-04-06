def rename_hits(species_name, out_dir, transcripts, db, basename):
    print('Renaming hits')
    length_seqs = [len(seq_record2.seq) for seq_record2 in
                   SeqIO.parse(db, "fasta")]  # create list with database sequences' lengths
    mean_length = statistics.mean(length_seqs)  # calculate mean of sequences' lengths
    my_records = []  # create list for writing of renamed hits
    filename = f'{out_dir}/{species_name}_hits.fasta'  # define file name of translated hits (merged with with initial filename)
    for seq_record in SeqIO.parse(f'{out_dir}/{basename}_blast_hits_clust.fasta', "fasta"):  # parse file with blast hits
        if transcripts == "cds":  # if cds mode has been chosen - rename hits and filter them
            if str(seq_record.seq)[0] == 'M' and len(str(seq_record.seq)) >= \
                    mean_length // 2:  # if seq starts with Met and its len is longer than 1/2 of mean len
                name = seq_record.id[:seq_record.id.find(" ")]  # take part of record name before space
                final = f'{species_name}_{name}'  # define final name of sequence
                rec = SeqRecord(seq_record.seq, id=final, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with sequence records
        elif transcripts == "all":  # if all mode has been chosen - just rename hits
            name = seq_record.id[:seq_record.id.find(" ")]
            print(name)
            final = f'{species_name}_{name}'  # define final name of sequence
            rec = SeqRecord(seq_record.seq, id=final, description='')
            my_records.append(rec)
    SeqIO.write(my_records, filename, 'fasta')
    print('Renaming is completed')
    return filename