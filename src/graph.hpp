#pragma once
#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <vector>
#include <climits>

struct Vertex {
  int n;         // The index of this vertex.
  //bool visited;  // Field for DFS/BFS.

  //int rank;  // field for treedepth on tree calculation

  //Vertex(int n) : n(n), visited(false) {}
  Vertex(int n) : n(n) {}
};

struct Graph {
  int N;                              // Number of vertices in this graph
  int M;                              // Number of edges in this graph.
  std::vector<Vertex> vertices;       // Vector containing all vertices.
  std::vector<std::vector<int>> adj;  // The adjacency lists for the full graph.

  Graph(std::istream &stream);
  Graph() : N(0), M(0) {}
};

struct Separator {
  std::vector<int> vertices;

  std::vector<std::pair<int, int>> comp;

  int maxCompSize() const {
    int result = 0;
    for (auto [N, M] : comp) result = std::max(result, N);
    return result;
  }
};

struct SubGraph {
  size_t max_degree = 0;        // Max degree of nodes inside this graph.
  size_t min_degree = INT_MAX;  // Min degree of nodes inside this graph.
  int N = 0;                    // Number of vertices in this graph.
  int M = 0;                    // Number of edges in this graph.

  std::vector<int> global;      // The global coordinates of the vertices in
                                // this graph.
  std::vector<std::vector<int>> adj;  // Adjacency list (local indexing).

  // Create an empty SubGraph.
  SubGraph();

  // Create a SubGraph of G with the given (local) vertices
  SubGraph(const SubGraph &G, const std::vector<int> &sub_vertices);

  // Checks whether this really represents an induced subgraph.
  // Note: expensive!
  void AssertValidSubGraph() const;

  // Vector of connected components of the subset given by sub_vertices.
  std::vector<SubGraph> ConnectedSubGraphs(
      const std::vector<int> &sub_vertices) const;

  // Get the adjacency list for a given vertex.
  const std::vector<int> &Adj(int v) const;

  // Get the local coordinate for a given vertex.
  int LocalIndex(int global_index) const;

  // Get all minimal separators for the given graph (as lists of local
  // coordinates).
  std::vector<Separator> AllMinimalSeparators() const;

  bool FullyMinimal(Separator &) const;

  // Create a connected components of the subgraph without the given vertices.
  std::vector<SubGraph> WithoutVertices(const std::vector<int> &S) const;

  // Create a connected components of the subgraph without the given vertex.
  std::vector<SubGraph> WithoutVertex(int v) const;

  // Recursively removes all vertices with deg < 2.
  SubGraph TwoCore() const;
  std::vector<SubGraph> kCore(int k) const;

  // Do a BFS from the given vertex.
  std::vector<int> Bfs(int v) const;

  // Compute trees from the given roots.
  SubGraph BfsTree(int root) const;
  SubGraph DfsTree(int root) const;

  // Returns whether this is a complete graph.
  bool IsCompleteGraph() const {
    return N * (N - 1) == 2 * M;
  }

  // Returns whether this is a path graph.
  bool IsPathGraph() const {
    return (N - 1 == M) && (max_degree < 3);
  }

  // Returns whether this is a star graph.
  bool IsStarGraph() const {
    return (N - 1 == M) && (M == max_degree);
  }

  // Returns whether this is a cycle graph.
  bool IsCycleGraph() const {
    return (M == N) && (max_degree == 2);
  }

  // Returns whether this is a tree.
  bool IsTreeGraph() const {
    return N - 1 == M;
  }

  // Explicit conversion to vector of ints.
  operator std::vector<int>() const {
    std::vector<int> result = global;
    std::sort(result.begin(), result.end());
    return result;
  }
};

extern Graph full_graph;  // The datastructure containing the full graph.
extern SubGraph full_graph_as_sub;  // The full graph in a SubGraph format.

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
