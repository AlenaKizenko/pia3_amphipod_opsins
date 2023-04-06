def match_amino_contig(filename, transcriptome, opsins = False):
    print('Matching protein sequence to contig in transcriptome')
    file_hits_set = {title for title, seq in SimpleFastaParser(
            open(f'{os.path.splitext(filename)[0]}_PIA3_aa.fasta'))}  # create set from seq names from classified opsins' file
    my_records = []  # create list for
    for seq_record in SeqIO.parse(transcriptome, "fasta"):  # parse transcriptome file
        for opsin in file_hits_set:  # iterate on names from classified opsins' file
            if seq_record.id.split()[0][:-1] in opsin:  # if seq names are equal
                rec = SeqRecord(seq_record.seq, id=opsin, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, f'{os.path.splitext(filename)[0]}_PIA3_contigs.fasta', 'fasta')  # write seq records to file
    print('Matching protein sequences to contigs is completed')