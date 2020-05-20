"""
usage: python3 vis.py outputfile.tree [inputfile.gr]
    
       Plots the tree specified in the outpufile, 
       if inputfile is specified it plots the tree but with edges from the original graph
"""

import igraph as ig
import sys

graphnr = sys.argv[1]


# reading treedepth decomposition
f = open("exact/exact_{}.tree".format(graphnr))
lines = f.readlines()

n = len(lines)-1
td = int(lines[0])

parents = [int(x)-1 for x in lines[1:]]
children = [[] for _ in range(n)]
edges = []
root = 0
for child, parent in enumerate(parents):
    if parent != -1:
        edges.append([parent, child])
        children[parent].append(child)
    else:
        root = child 

# determine first separator, for sure the root is included, 
# then keep on adding nodes while the current node only has a single child
first_separator = [root] 
if len(children[root]) == 1:
    cur_node = root
    while len(children[cur_node]) == 1:
        cur_node = children[cur_node][0]
        first_separator.append(cur_node)

# vertex colors (green if in the first separator otherwise red 
vertex_color = ["green" if x in first_separator else "red" for x in range(n)]
# vertex labels (index + 1)
vertex_labels = [str(x+1) for x in range(n)]

g = ig.Graph(n = n, directed = True)
g.add_edges(edges)

g.vs["label"] = vertex_labels 
g.vs["color"] = vertex_color

treelayout = g.layout_reingold_tilford(root=root)

f.close()

# reading input graph
g2 = ig.Graph(n=n)
f = open("../input/exact/exact_{}.gr".format(graphnr))
line = f.readline().split()
m = int(line[3])
for _ in range(m):
    u, v = [int(x) for x in f.readline().split()]
    g2.add_edge(u-1, v-1)

g2.vs["label"] = vertex_labels 
g2.vs["color"] = vertex_color


ig.plot(g, layout = treelayout)
#ig.plot(g2, layout = treelayout)

#layout2 = g2.layout("drl")
#layout2 = g2.layout("grid_fr")
layout2 = g2.layout("large")
ig.plot(g2, layout = layout2)





