#pragma once
#include <iostream>
#include <vector>

struct Vertex {
  int n_;                      // The index of this vertex.
  std::vector<Vertex *> adj_;  // The adjacency list of this vertex.
  bool visited_;               // Field for DFS/BFS.

  Vertex(int n) : n_(n) {}
};

struct Graph {
  int N_;                         // Number of vertices in this graph
  int M_;                         // Number of edges in this graph.
  std::vector<Vertex> vertices_;  // Vector containing all vertices.

  Graph(std::istream &stream);
  Graph() : N_(0), M_(0) {}
};

struct SubGraph {
  std::vector<Vertex *> vertices_;  // List of vertices inside this subgraph.
  std::vector<bool> mask_;  // Bitset of the vertices inside this subgraph.

  SubGraph(std::vector<Vertex *> &&vertices, std::vector<bool> &&mask);
  SubGraph();

  // Create a new subgraph withouth the given vertex.
  SubGraph WithoutVertex(Vertex *v);
  std::vector<SubGraph> ConnectedComponents();
};

extern Graph full_graph;  // The datastructure containing the full graph.
extern SubGraph full_graph_as_sub;  // The full graph in a SubGraph format.

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
