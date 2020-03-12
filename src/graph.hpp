#pragma once
#include <algorithm>
#include <climits>
#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <vector>

struct Vertex {
  int n;         // The index of this vertex.
  bool visited;  // Field for DFS/BFS.

  int rank;  // field for treedepth on tree calculation

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
  size_t max_degree = 0;        // Max degree of nodes inside this graph.
  size_t min_degree = INT_MAX;  // Min degree of nodes inside this graph.
  int M = 0;                    // Number of edges in this subgraph.

  std::vector<Vertex *> vertices;  // List of vertices inside this subgraph.
  std::vector<bool> mask;  // Bitset of the vertices inside this subgraph.
  std::vector<std::vector<int>> adj;  // Adjacency list in local indexing.

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
  int LocalIndex(Vertex *v) const;

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

  // Returns whether this is a cycle graph.
  bool IsCycleGraph() const {
    int N = vertices.size();
    return (M == N) && (max_degree == 2);
  }

  // Returns whether this is a tree.
  bool IsTreeGraph() const {
    int N = vertices.size();
    return N - 1 == M;
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

struct Separator {
  std::vector<int> vertices;
  std::vector<std::pair<int, int>> comp;

  Separator() {}
  Separator(std::vector<int> &&vertices) : vertices(std::move(vertices)) {}

  int maxCompSize() const {
    int result = 0;
    for (auto [N, M] : comp) result = std::max(result, N);
    return result;
  }
};

class SeparatorGenerator {
 public:
  SeparatorGenerator(const SubGraph &G);

  bool HasNext() const { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000);

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // Helper function.
  bool FullyMinimal(Separator &) const;

  const SubGraph &G;
  int N;

  // In done we keep the seperators we have already enqueued, to make sure they
  // aren't processed again.
  // In queue we keep all the ones we have generated, but which we have not yet
  // used to generate new ones.
  std::queue<std::vector<int>> queue;
  std::set<std::vector<int>> done;

  std::vector<bool> in_nbh;
};

extern Graph full_graph;  // The datastructure containing the full graph.
extern SubGraph full_graph_as_sub;  // The full graph in a SubGraph format.

// This initalizes the above global variables, important!
void LoadGraph(std::istream &stream);
