#!/usr/bin/env python3
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
T = ['DendroPy','Bio.Phylo','TreeSwift']
pal = {'TreeSwift':'#0000FF', 'Bio.Phylo':'#00FF00', 'DendroPy':'#FF0000'}
handles = [Patch(color=pal[t],label=t) for t in T]
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
    for t in T:
        x = []; y = []
        for n in N:
            if n not in data[m]:
                continue
            if t.lower().replace('.','') in data[m][n]:
                x += [n]*10
                y += data[m][n][t.lower().replace('.','')]
            else:
                continue
        if len(x) != 0:
            x = np.array(x); y = np.array(y)
            sns.pointplot(x=x, y=y, color=pal[t])
    plt.ylim(ymin=0)
    sns.plt.title("Execution Time vs. Number of Leaves (%s)"%task[m],fontsize=18,y=1.02)
    sns.plt.xlabel("Number of Leaves")
    sns.plt.ylabel("Execution Time (seconds)")
    legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
    plt.show()
    fig.savefig('%s.pdf'%m.lower(), format='pdf', bbox_inches='tight')
