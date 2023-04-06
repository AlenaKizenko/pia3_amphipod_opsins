def classify_opsins(tree, species):
    print('Classifying opsins by sensitivity')
    global opsins_class  # assign global variable
    opsins_class = {}  # create dictionary for opsins' types classification
    for leaf in tree.iter_leaves():  # iterate on leaves
        if 'RHO' in str(leaf):  # if RHO (Outgroup - OG) is in reference database sequence name
            while leaf.up:  # move to the root of the tree
                if not str(leaf.name):  # if node name is empty
                    leaf.name = 'OG'  # name node as OG
                leaf = leaf.up  # move to the root of the tree once
        elif 'MEL' in str(leaf) or 'PER' in str(leaf) or 'TMT' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'VERL'
                leaf = leaf.up
        elif 'MWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'MWS'
                leaf = leaf.up
        elif 'LWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'LWS'
                leaf = leaf.up
        elif 'SWS' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'SWS'
                leaf = leaf.up
        elif 'UV' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'UV'
                leaf = leaf.up

    for leaf in tree.iter_leaves():  # iterate on leaves
        if species in str(leaf):  # if leaf is a hit
            parent = leaf.up  # parent is a node
            while not str(parent.name):  # while node name is empty
                parent = parent.up  # move to the root of the tree
            type_opsin = str(parent.name)  # define type of opsin as parent node name
            opsins_class[leaf.name] = type_opsin  # write type of hit to dict
            leaf.name = f'{type_opsin}_{leaf.name}'  # rename leaf name
    print('Classification is completed')
    return tree