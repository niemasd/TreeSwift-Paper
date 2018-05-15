#!/usr/bin/env python3
# get memory usage
def memory():
    from os import getpid; from psutil import Process
    return Process(getpid()).memory_info().rss

# main code
from sys import argv
from time import time
if len(argv) != 3 or argv[2] not in {'treeswift','dendropy'}:
    print("USAGE: %s <tree_file> <treeswift_or_dendropy>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()
m_start = memory(); t_start = time()
if argv[2][0] == 'd':
    from dendropy import Tree
    for node in Tree.get(data=treestr, schema='newick').preorder_node_iter():
        pass
else:
    from treeswift import read_tree_newick
    for node in read_tree_newick(treestr).traverse_preorder():
        pass
t_end = time(); m_end = memory()
print("Time: %f"%(t_end-t_start))
print("Memory: %d"%(m_end-m_start))
