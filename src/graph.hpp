#pragma once
#include <algorithm>
#include <cassert>
#include <climits>
#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <vector>

struct Graph {
  size_t max_degree = 0;        // Max degree of nodes inside this graph.
  size_t min_degree = INT_MAX;  // Min degree of nodes inside this graph.
  int N = 0;                    // Number of vertices in this graph.
  int M = 0;                    // Number of edges in this graph.

  std::vector<int> global;  // The global coordinates of the vertices in
                            // this graph.
  std::vector<std::vector<int>> adj;  // Adjacency list (local indexing).

  // Create an empty Graph.
  Graph();

  // Create a Graph from a stream.
  Graph(std::istream &stream);

  // Create a Graph of G with the given (local) vertices.
  Graph(const Graph &G, std::vector<int> &sub_vertices);

  // Checks whether this really represents an induced subgraph.
  // Note: expensive!
  void AssertValidGraph() const;

  // Vector of connected components of the subset given by sub_vertices.
  std::vector<Graph> ConnectedGraphs(
      const std::vector<int> &sub_vertices) const;

  // Get the adjacency list for a given vertex.
  inline const std::vector<int> &Adj(int v) const {
    assert(v >= 0 && v < N && adj.size() == N);
    return adj[v];
  }

  // Get the local coordinate for a given vertex.
  inline int LocalIndex(int global_index) const {
    for (int v_local = 0; v_local < N; ++v_local) {
      if (global[v_local] == global_index) return v_local;
    }
    assert(false);
  }

  // Checks if the subset in vertices is a connected subset of the graph.
  bool ConnectedSubset(const std::vector<int> vertices) const;

  // Contract the vertices in `contractors` into a single vertex. Assumes that
  // `contractors` forms a connected subset of vertices.
  Graph Contract(const std::vector<int> &contractors) const;

  // Create a connected components of the subgraph without the given vertices.
  std::vector<Graph> WithoutVertices(const std::vector<int> &S) const;

  // Create a connected components of the subgraph without the given vertex.
  std::vector<Graph> WithoutVertex(int v) const;

  // Removes all vertices w > v for which  N(w) = N(v) or N(w)\v = N(v)\w.
  // It returns the smaller, and for each vertex in the new graph
  // a list of vertices in the old graph that share the same neighboorhood.
  std::pair<Graph, std::vector<std::vector<int>>>
  WithoutSymmetricNeighboorhoods() const;

  // Returns a sublist of the vertices that are not dominated by a vertex
  // in the list.
  std::vector<int> NonDominatedVertices(const std::vector<int> &vertices) const;

  // Recursively removes all vertices with deg < 2.
  Graph TwoCore() const;
  std::vector<Graph> kCore(int k) const;

  // Do a BFS from the given vertex.
  std::vector<int> Bfs(int v) const;

  // Compute trees from the given roots.
  Graph BfsTree(int root) const;
  Graph DfsTree(int root) const;

  // Computes a list of all articulation points.
  std::vector<int> ArticulationPoints() const;

  // Returns whether this is a complete graph.
  inline bool IsCompleteGraph() const { return N * (N - 1) == 2 * M; }

  // Returns whether this is a path graph.
  inline bool IsPathGraph() const { return (N - 1 == M) && (max_degree < 3); }

  // Returns whether this is a star graph.
  inline bool IsStarGraph() const { return (N - 1 == M) && (M == max_degree); }

  // Returns whether this is a cycle graph.
  inline bool IsCycleGraph() const { return (M == N) && (max_degree == 2); }

  // Returns whether this is a tree.
  inline bool IsTreeGraph() const { return N - 1 == M; }

  // Explicit conversion to vector of ints.
  operator std::vector<int>() const {
    std::vector<int> result = global;
    std::sort(result.begin(), result.end());
    return result;
  }
};

extern Graph full_graph;                   // The full graph.
extern std::vector<bool> full_graph_mask;  // Global variable to be reused.

// For going from global coordinates to sets of original vertices, and back.
extern std::vector<std::vector<int>> global_to_vertices;
extern std::map<std::vector<int>, int> vertices_to_global;

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
