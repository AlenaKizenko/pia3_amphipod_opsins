def check_lysine(out_dir, alignment, filename, ref_seq_name = 'RHO_Bos_taurus_AAA30674.1', n=296): # calculate position in the alignment
    print('Checking lysine position in opsins\' sequences')
    with open(alignment) as aln:
        for record in AlignIO.read(aln, "fasta"):
            if str(record.id) == ref_seq_name:  # if the ref sequence found
                sequence_str = record.seq
                letter_counter = 0
                gap_counter = 0
                for letter in sequence_str:
                    if letter_counter < n - 1:
                        if letter == "-":
                            gap_counter += 1
                        else:
                            letter_counter += 1
                break

    number_position = gap_counter + letter_counter

    all_opsins = []
    with open(alignment) as aln:
        for seq_record in AlignIO.read(alignment, "fasta"):
            if str(seq_record.seq)[number_position] == 'K':
                all_opsins.append(seq_record.id)       
    query_opsins = []
    for seq_record in SeqIO.parse(f'{out_dir}/{filename}_hits.fasta', 'fasta'):
        for opsin in all_opsins:
            if opsin in seq_record.id:
                query_opsins.append(seq_record)
    SeqIO.write(query_opsins, f'{filename}_opsins.fasta', 'fasta')
    print('Checking lysine position is completed')
    return f'{filename}_opsins.fasta'