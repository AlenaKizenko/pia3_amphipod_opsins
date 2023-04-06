hits = set()  # create set for blast hits: since names can be repeated, there is no need to use list
    my_records = []  # create list for hits to write them to .fasta file
    with open(f'{out_dir}/{basename}_blast_file.tmp') as blast_hits:  # open blast file
        for blast_hit in blast_hits:
            hits.add(blast_hit.split()[
                         0])  # write blast hits' name to the file (split because of Transdecoder's prolonged names)
    for seq_record in SeqIO.parse(cds, "fasta"):
        for hit in hits:  # iterate on all transcriptes and blast hits
            if str(hit) == str(seq_record.id):  # if names are equal
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # make SeqRecord object
                my_records.append(rec)  # write seq record to list my_records
    SeqIO.write(my_records, f'{out_dir}/{basename}_blast_hits.fasta', 'fasta')  # write records to .fasta file
    print('BLAST search is completed')