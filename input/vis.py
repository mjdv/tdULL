"""
usage: python3 vis.py inputfile.gr. Plots the graph specified in the inputfile. 
"""

import sys

import igraph as ig

f = open(sys.argv[1])
line = f.readline().split()
n = int(line[2])
m = int(line[3])

g = ig.Graph(n=n, directed=False)
for _ in range(m):
    u, v = [int(x) for x in f.readline().split()]
    g.add_edge(u - 1, v - 1)
g.vs["label"] = [str(x) for x in range(n)]
ig.plot(g)
