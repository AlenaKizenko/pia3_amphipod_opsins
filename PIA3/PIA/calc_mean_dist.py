def calc_mean_dist(tree):
    print('Calculating median absolute deviation of evolutionary distances')
    tree = ete3.Tree(tree)  # assign Ete3 tree object
    lst = []  # create list for distances
    for leaf in tree.iter_leaves():  # iterate on leaves
        dist = leaf.dist  # assign variable for distance
        lst.append(dist)  # append list with distances
    me = statistics.mean(lst)  # calculate distances mean
    lst_me = []  # create list for absolute deviations from mean
    for i in lst:
        a = abs(i - me)  # from each distance subtract mean and take absolute value
        lst_me.append(a)  # append list with absolute deviation from mean
    print('Calculations are completed')
    return statistics.mean(lst_me) * 4  # return mean absolute deviation * 4