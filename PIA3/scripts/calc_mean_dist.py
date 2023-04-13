import ete3
import statistics
import sys

def calc_mean_dist(tree):
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
    return statistics.mean(lst_me) * 4  # return mean absolute deviation * 4

tree = sys.argv[1]
output_file = sys.argv[2]

calc_result = calc_mean_dist(tree)

with open(output_file, 'w') as out:
    out.write(f'{calc_result}')

