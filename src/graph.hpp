#pragma once
#include <iostream>
#include <vector>

struct Vertex {
  int n;                      // The index of this vertex.
  std::vector<Vertex *> adj;  // The adjacency list of this vertex.
  bool visited;               // Field for DFS/BFS.

  Vertex(int n) : n(n) {}
};

struct Graph {
  int N;                         // Number of vertices in this graph
  int M;                         // Number of edges in this graph.
  std::vector<Vertex> vertices;  // Vector containing all vertices.

  Graph(std::istream &stream);
  Graph() : N(0), M(0) {}
};

struct SubGraph {
  std::vector<Vertex *> vertices;  // List of vertices inside this subgraph.
  std::vector<bool> mask;  // Bitset of the vertices inside this subgraph.

  SubGraph(std::vector<Vertex *> &&vertices, std::vector<bool> &&mask);
  SubGraph();

  // Get the adjacency for a vertex inside this subgraph.
  std::vector<Vertex *> Adj(Vertex *v);

  // Create a new subgraph withouth the given vertex.
  SubGraph WithoutVertex(Vertex *v);

  // Get all the connected components.
  std::vector<SubGraph> ConnectedComponents();
};

extern Graph full_graph;  // The datastructure containing the full graph.
extern SubGraph full_graph_as_sub;  // The full graph in a SubGraph format.

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
