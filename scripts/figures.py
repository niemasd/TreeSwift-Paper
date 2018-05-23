#!/usr/bin/env python3
MEMORY_DEN = 1024*1024

# load data
from pickle import load
from sys import argv
if len(argv) != 2:
    print("USAGE: %s <data.pkl>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    data = load(gopen(argv[1]))
else:
    data = load(open(argv[1]))

# set up seaborn
from matplotlib import rcParams
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style("ticks")
rcParams['font.family'] = 'serif'
T = ['DendroPy','Bio.Phylo','ETE Toolkit','TreeSwift']
name_to_key = {'DendroPy':'dendropy', 'Bio.Phylo':'biophylo', 'TreeSwift':'treeswift', 'ETE Toolkit':'ete3'}
pal = {'TreeSwift':'#20641c', 'Bio.Phylo':'#6caed1', 'DendroPy':'#2286ca', 'ETE Toolkit':'#bfe49e'}
linestyle = {'TreeSwift':'-', 'DendroPy':':', 'Bio.Phylo':'--', 'ETE Toolkit':'-.'}
task = {
    'distance_matrix': "Distance Matrix",
    'inorder': "In-Order",
    'levelorder': "Level-Order",
    'mrca': "MRCA",
    'postorder': "Post-Order",
    'preorder': "Pre-Order",
    'rootdistorder': "Root-Distance-Order",
    "total_branch_length": "Total Length"
}
N = None
for m in data:
    if N is None or len(data[m]) > len(N):
        N = sorted(data[m].keys())

# create figures
for m in sorted(data.keys()):
    fig = plt.figure()
    used = []
    for t in T:
        x = []; y = []
        for n in N:
            if n not in data[m]:
                continue
            if name_to_key[t] in data[m][n]:
                x += [n]*10
                y += data[m][n][name_to_key[t]]
            else:
                continue
        if len(x) != 0:
            x = np.array(x); y = np.array(y)
            if m == 'memory':
                y = y/float(MEMORY_DEN)
            sns.pointplot(x=x, y=y, color=pal[t], linestyles=linestyle[t])
            used.append(t)
    plt.ylim(ymin=0)
    if m == 'memory':
        sns.plt.title("Memory Usage vs. Number of Leaves",fontsize=18,y=1.02)
    elif m == 'rootdistorder':
        sns.plt.title("Execution Time vs. Number of Leaves (%s)"%task[m],fontsize=14,y=1.02)
    else:
        sns.plt.title("Execution Time vs. Number of Leaves (%s)"%task[m],fontsize=16,y=1.02)
    sns.plt.xlabel("Number of Leaves")
    if m == 'memory':
        sns.plt.ylabel("Memory Usage (MB)")
    else:
        sns.plt.ylabel("Execution Time (seconds)")
    handles = [Patch(color=pal[t],label=t) for t in used]
    legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
    plt.show()
    fig.savefig('%s.pdf'%m.lower(), format='pdf', bbox_inches='tight')