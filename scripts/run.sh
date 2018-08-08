#!/usr/bin/env bash
IT=10
m=$1
n=$2
t=$3
echo -n "["
for r in $(seq -w 1 $IT) ; do ./time.py "../data/tree_n$n.tre.gz" $t $m ; done | numlist -csv | tr -d '\n'
echo "]"
