#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

struct Vertex {
  int n;         // The index of this vertex.
  bool visited;  // Field for DFS/BFS.

  int rank; // field for treedepth on tree calculation

  Vertex(int n) : n(n), visited(false) {}
};

struct Graph {
  int N;                              // Number of vertices in this graph
  int M;                              // Number of edges in this graph.
  std::vector<Vertex> vertices;       // Vector containing all vertices.
  std::vector<std::vector<int>> adj;  // The adjacency lists for the full graph.

  Graph(std::istream &stream);
  Graph() : N(0), M(0) {}
};

struct SubGraph {
  size_t max_degree = 0;  // Max degree of nodes inside this graph.
  int M = 0;              // Number of edges in this subgraph.

  std::vector<Vertex *> vertices;  // List of vertices inside this subgraph.
  std::vector<bool> mask;  // Bitset of the vertices inside this subgraph.
  std::vector<std::vector<int>> adj;  // Adjacency list in local indexing.

  // Create an empty SubGraph.
  SubGraph();

  // Get the adjacency list for a given vertex.
  const std::vector<int> &Adj(int v) const;

  // Create a connected components of the subgraph without the given vertex.
  std::vector<SubGraph> WithoutVertex(int v) const;

  // Do a BFS from the given vertex.
  std::vector<int> Bfs(int v) const;

  // Returns whether this is a complete graph.
  bool IsCompleteGraph() const {
    int N = vertices.size();
    return N * (N - 1) == 2 * M;
  }

  // Returns whether this is a path graph.
  bool IsPathGraph() const {
    int N = vertices.size();
    return (N - 1 == M) && (max_degree < 3);
  }

  // Returns whether this is a star graph.
  bool IsStarGraph() const {
    int N = vertices.size();
    return (N - 1 == M) && (M == max_degree);
  }

  // Explicit conversion to vector of ints.
  operator std::vector<int>() const {
    std::vector<int> result;
    result.reserve(vertices.size());
    for (Vertex *v : vertices) result.emplace_back(v->n);
    std::sort(result.begin(), result.end());
    return result;
  }
};

extern Graph full_graph;  // The datastructure containing the full graph.
extern SubGraph full_graph_as_sub;  // The full graph in a SubGraph format.

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
