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
N = [100,1000,10000,100000,1000000]
T = ['TreeSwift','Bio.Phylo','DendroPy']
units = {'Time':'s','Memory':'MB'}
pal = {'TreeSwift':'#0000FF', 'Bio.Phylo':'#00FF00', 'DendroPy':'#FF0000'}
handles = [Patch(color=pal[t],label=t) for t in T]

# create figures
for m in ['Time','Memory']:
    fig = plt.figure()
    for t in T:
        x = []; y = []
        for n in N:
            x += [n]*10; y += data[t.lower().replace('.','')][m.lower()][n]
        x = np.array(x); y = np.array(y)
        if m == 'Memory':
            y /= (1024*1024)
        sns.pointplot(x=x, y=y, color=pal[t])
        sns.plt.title("%s vs. Number of Leaves"%m,fontsize=18,y=1.02)
        sns.plt.xlabel("Number of Leaves")
        sns.plt.ylabel("%s (%s)"%(m,units[m]))
        legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
    plt.show()
    fig.savefig('%s_vs_n.pdf'%m.lower(), format='pdf', bbox_inches='tight')
