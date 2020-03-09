"""
usage: python3 vis.py outputfile.tree [inputfile.gr]
    
       Plots the tree specified in the outpufile, 
       if inputfile is specified it plots the tree but with edges from the original graph
"""

import igraph as ig
import sys

f = open(sys.argv[1])
lines = f.readlines()

n = len(lines)-1
td = int(lines[0])

parents = [int(x) for x in lines[1:]]
edges = []
root = 0
for i, p in enumerate(parents):
    if p != 0:
        edges.append([p-1, i])
    else:
        root = i


g = ig.Graph(n = n, directed = True)
g.add_edges(edges)

g.vs["label"] = [str(x+1) for x in range(n)]

layout = g.layout_reingold_tilford(root=root)

f.close()

if len(sys.argv) == 3:
    g.delete_edges(ig.EdgeSeq(g))
    f = open(sys.argv[2])
    line = f.readline().split()
    m = int(line[3])
    for _ in range(m):
        u, v = [int(x) for x in f.readline().split()]
        g.add_edge(u-1, v-1)

    g.to_undirected()

ig.plot(g, layout = layout)







