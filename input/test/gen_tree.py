import sys
import networkx as nx

n = int(sys.argv[1])

g = nx.random_tree(n)

f = open("{}_random_tree.gr".format(n), "w")
f.write("p tdp {} {}\n".format(n, n-1))
for u,v in g.edges():
    f.write("{} {}\n".format(u+1, v+1))

f.close()
