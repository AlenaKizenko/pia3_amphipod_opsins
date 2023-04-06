def filter_distant_seqs(input_file, tree_query, dist_dev, query_file):
    print('Filtering distinct hits')
    tree_query = ete3.Tree(tree_query)  # assign Ete3 tree object
    lst_seqs = []  # create list for selected leaves (distance less than 4*absolute mean deviation)
    my_records = []  # create list for selected sequences
    for leaf in tree_query.iter_leaves():  # iterate on leaves of query tree
        if input_file in str(
                leaf.name) and leaf.dist < dist_dev:  # if seq name from file == seq name from tree and distance is OK
            lst_seqs.append(str(leaf.name))  # append lst_seqs by name of selected leaf
    for seq_record in SeqIO.parse(query_file, "fasta"):  # parse .fasta file with hits
        for seq_name in lst_seqs:
            seq_record.id = seq_record.id.replace(":", "_")  # replacement because of IQ-Tree specificity
            seq_record.id = seq_record.id.replace("|", "_")  # replacement because of IQ-Tree specificity
            if seq_record.id in seq_name:  # if seq name from file == selected name
                rec = SeqRecord(seq_record.seq, seq_record.id, description='')  # create SeqRecord object
                my_records.append(rec)  # append list with seq records
    SeqIO.write(my_records, f'{input_file}_PIA3_aa.fasta', 'fasta')  # write seq records to file
    print('Filtering is completed')