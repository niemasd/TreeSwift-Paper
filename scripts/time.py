#!/usr/bin/env python3
from Bio import Phylo
from dendropy import Tree
from treeswift import read_tree_newick
import numpy

# main code
from io import StringIO
from sys import argv
from time import time
if len(argv) != 4 or argv[2] not in {'treeswift','dendropy','biophylo'}:
    print("USAGE: %s <tree_file> <treeswift_or_dendropy_or_biophylo> <task>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()
treeio = StringIO(treestr) # for Bio.Phylo

# distance matrix
def distance_matrix(m):
    if m == 'dendropy':
        t_start = time()
        Tree.get(data=treestr, schema='newick').phylogenetic_distance_matrix()
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        # from http://biopython.org/wiki/Phylo_cookbook
        tree = Phylo.read(treeio, 'newick')
        allclades = list(tree.find_clades(order='level'))
        lookup = {}
        for i, elem in enumerate(allclades):
            lookup[elem] = i
        distmat = numpy.repeat(numpy.inf, len(allclades)**2)
        distmat.shape = (len(allclades), len(allclades))
        for parent in tree.find_clades(terminal=False, order='level'):
            for child in parent.clades:
                if child.branch_length:
                    distmat[lookup[parent], lookup[child]] = child.branch_length
        if not tree.rooted:
            distmat = distmat + distmat.transpose()
        numpy.matrix(distmat)
        t_end = time()
    else:
        t_start = time()
        read_tree_newick(treestr).distance_matrix()
        t_end = time()
    return t_end-t_start

# inorder traversal
def inorder(m):
    if m == 'dendropy':
        t_start = time()
        for node in Tree.get(data=treestr, schema='newick').inorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        raise RuntimeError("Bio.Phylo does not have this function")
    else:
        t_start = time()
        for node in read_tree_newick(treestr).traverse_inorder():
            pass
        t_end = time()
    return t_end-t_start

# level-order traversal
def levelorder(m):
    if m == 'dendropy':
        t_start = time()
        for node in Tree.get(data=treestr, schema='newick').levelorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        for node in Phylo.read(treeio, 'newick').find_clades(order='level'):
            pass
        t_end = time()
    else:
        t_start = time()
        for node in read_tree_newick(treestr).traverse_levelorder():
            pass
        t_end = time()
    return t_end-t_start

# postorder traversal
def postorder(m):
    if m == 'dendropy':
        t_start = time()
        for node in Tree.get(data=treestr, schema='newick').postorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        for node in Phylo.read(treeio, 'newick').find_clades(order='postorder'):
            pass
        t_end = time()
    else:
        t_start = time()
        for node in read_tree_newick(treestr).traverse_postorder():
            pass
        t_end = time()
    return t_end-t_start

# preorder traversal
def preorder(m):
    if m == 'dendropy':
        t_start = time()
        for node in Tree.get(data=treestr, schema='newick').preorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        for node in Phylo.read(treeio, 'newick').find_clades(order='preorder'):
            pass
        t_end = time()
    else:
        t_start = time()
        for node in read_tree_newick(treestr).traverse_preorder():
            pass
        t_end = time()
    return t_end-t_start

# root distance order traversal
def rootdistorder(m):
    if m == 'dendropy':
        t_start = time()
        tree = Tree.get(data=treestr, schema='newick')
        tree.calc_node_ages(is_force_max_age=True)
        for node in tree.ageorder_node_iter(descending=True):
            pass
        t_end = time()
    elif m == 'biophylo':
        raise RuntimeError("Bio.Phylo does not have this function")
    else:
        t_start = time()
        for node in read_tree_newick(treestr).traverse_rootdistorder():
            pass
        t_end = time()
    return t_end-t_start

# total branch length
def total_branch_length(m):
    if m == 'dendropy':
        t_start = time()
        Tree.get(data=treestr, schema='newick').length()
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        Phylo.read(treeio, 'newick').total_branch_length()
        t_end = time()
    else:
        t_start = time()
        read_tree_newick(treestr).edge_length_sum()
        t_end = time()
    return t_end-t_start

TASKS = {
    'distance_matrix':distance_matrix,
    'inorder':inorder,
    'levelorder':levelorder,
    'postorder':postorder,
    'preorder':preorder,
    'rootdistorder':rootdistorder,
    'total_branch_length':total_branch_length,
}

# run
if argv[3] not in TASKS:
    print("Invalid task. Valid options: %s"%', '.join(sorted(TASKS.keys()))); exit(1)
print(TASKS[argv[3]](argv[2]))