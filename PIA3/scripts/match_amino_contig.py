def match_amino_contig(input_file, transcriptome, output_file):
    file_hits_set = {title for title, seq in SimpleFastaParser(
            open(input_file)}  # create set from seq names from classified opsins' file
    my_records = []  # create list for
    for seq_record in SeqIO.parse(transcriptome, "fasta"):  # parse transcriptome file
        for hit in file_hits_set:  # iterate on names from classified opsins' file
            if seq_record.id.split()[0][:-1] in hit:  # if seq names are equal
                rec = SeqRecord(seq_record.seq, id=hit, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, output_file, 'fasta')  # write seq records to file


input_file = sys.argv[1]
transcriptome = sys.argv[2]
output_file = sys.argv[3]


result_cds = match_amino_contig(input_file, transcriptome, output_file)