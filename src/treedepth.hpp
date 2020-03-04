#pragma once
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>

#include "graph.hpp"
#include "set_trie.hpp"

// The global cache (a SetTrie) is what we use to store bounds on treedepths
// for subsets of the global graph (which correspond to induced subgraphs).
// Every subgraph that is in the cache comes with three data: a lower_bound, an
// upper_bound, and a root.
//
// The following is expect to be true of the cache and must be maintained:
// - It contains only connected subgraphs of the connected graph.
// - The lower_bound contains a proven lower bound on the subgraph.
// - The upper_bound contains a proven upper bound on the subgraph, where
// - The root is an element of the subgraph which witnesses this upper_bound,
//   and furthermore each connected component of the subgraph with the root
//   removed is also in the cache.
SetTrie cache;

// The function treedepth computes Treedepth bounds on subgraphs of the global
// graph.
//
// Parameters:
// - SubGraph G, the graph for which we try to compute treedepth bounds.
// - int search_lbnd, the lower bound of which treedepths are useful.
// - int search_ubnd, the upper bound of which treedepths are useful.
//
// Returns (as pair<int, int>):
// - lower: a lower bound on the treedepth of the graph.
// - upper: an upper bound on the treedepth of the graph.
//
// Explanation of the search upper bound: the instance of treedepth that has
// called this instance may already have a decomposition of depth d of the
// graph it is decomposing. Then there is no reason to continue the search for
// this subgraph if its treedepth lower bound exceeds d - 1 (= search_lbnd).
// Thus we can exit as soon as lower comes above search_ubnd.
//
// Explanation of the search lower bound: the instance of treedepth that has
// called this instance will possibly be calling it on multiple components. If
// one sister component already has a lower bound on its treedepth of d, then
// there is no reason to try to get the treedepth of this subgraph any lower
// than d. Thus if we find a decomposition that can has depth at most d, i.e.
// upper is at most search_lbnd, we are done.
std::pair<int, int> treedepth(const SubGraph &G, int search_lbnd,
                              int search_ubnd) {
  int N = G.vertices.size();
  assert(N >= 1);

  int lower = G.M / N + 1, upper = N;

  // Check whether this graph is in the cache.
  auto node = cache.Search(G);
  if (node != nullptr) {
    lower = node->data.lower_bound;
    upper = node->data.upper_bound;
  } else {
    node = cache.Insert(G);
    node->data.lower_bound = lower;
    node->data.upper_bound = upper;
    node->data.root = G.vertices[0]->n;
  }

  // If the trivial or previously found bounds suffice, we are done.
  if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
    return {lower, upper};
  }

  // Little helper function to add lower/upper/root to cache.
  auto CacheInsert = [node](int lower_bound, int upper_bound, int root) {
    node->data.lower_bound = lower_bound;
    node->data.upper_bound = upper_bound;
    node->data.root = root;
    return std::pair{lower_bound, upper_bound};
  };

  // We do a quick check for special cases we can answer exactly and
  // immediately.
  if (G.IsCompleteGraph()) return CacheInsert(N, N, G.vertices[0]->n);
  if (G.IsStarGraph()) {
    // Find node with max_degree.
    for (int v = 0; v < G.vertices.size(); ++v)
      if (G.Adj(v).size() == G.max_degree)
        return CacheInsert(2, 2, G.vertices[v]->n);
  }
  if (G.IsPathGraph()) {
    // Find the bound
    int bnd = 1;
    while (N >>= 1) bnd++;

    // Find a leaf, do BFS and find middle node.
    for (int v = 0; v < G.vertices.size(); ++v)
      if (G.Adj(v).size() == 1) {
        auto bfs = G.Bfs(v);
        return CacheInsert(bnd, bnd, G.vertices[bfs[bfs.size() / 2]]->n);
      }
  }

  // Create vector with numbers 0 .. N - 1
  std::vector<int> sorted_vertices(G.vertices.size());
  std::iota(sorted_vertices.begin(), sorted_vertices.end(), 0);

  // Sort the vertices based on the degree.
  std::sort(
      sorted_vertices.begin(), sorted_vertices.end(),
      [&](int v1, int v2) { return G.Adj(v1).size() > G.Adj(v2).size(); });

  // Main loop: try every vertex as root.
  // new_lower tries to find a new treedepth lower bound on this subgraph.
  int new_lower = N;

  // If the graph has at least 3 vertices, we never want a leaf (degree 1
  // node) as a root.
  bool skip_leaves = G.vertices.size() > 2;
  for (auto v : sorted_vertices) {
    if (skip_leaves && G.Adj(v).size() == 1) continue;
    int search_ubnd_v = std::min(search_ubnd - 1, upper - 1);
    int search_lbnd_v = std::max(search_lbnd - 1, 1);

    int upper_v = 0;
    int lower_v = lower - 1;

    bool early_break = false;

    for (auto H : G.WithoutVertex(v)) {
      auto [lower_H, upper_H] = treedepth(H, search_lbnd_v, search_ubnd_v);

      upper_v = std::max(upper_v, upper_H);
      lower_v = std::max(lower_v, lower_H);

      search_lbnd_v = std::max(search_lbnd_v, lower_H);

      if (lower_H >= search_ubnd_v) {
        // This component already shows that there's no reason to
        // continue trying with vertex v.
        early_break = true;
        break;
      }
    }

    new_lower = std::min(new_lower, lower_v + 1);

    // The upper bound we found for v is only meaningful if we didn't break
    // early.
    if (!early_break && upper_v + 1 < upper) {
      upper = upper_v + 1;
      node->data.upper_bound = upper;
      node->data.root = G.vertices[v]->n;
    }

    if (upper <= search_lbnd || lower == upper) {
      // Choosing root v already gives us a treedepth decomposition which is
      // good enough (either a sister branch is at least this long, or it
      // matches a previously proved lower bound for this subgraph) so we can
      // use v as our root.
      return {lower, upper};
    }
  }

  lower = std::max(lower, new_lower);
  node->data.lower_bound = lower;
  return {lower, upper};
}

// Recursive function to reconstruct the tree that atains the treedepth.
void reconstruct(const SubGraph &G, int root, std::vector<int> &tree) {
  assert(G.vertices.size());
  auto node = cache.Search(G);

  // Not all subgraphs are neccessarily inside the cache.
  if (node == nullptr) {
    treedepth(G, 1, G.vertices.size());
    node = cache.Search(G);
  }
  assert(node);
  assert(node->data.root > -1);
  tree.at(node->data.root) = root;

  // Root is the global coordinate, find its local coordinate.
  int local_root = -1;
  for (int v = 0; v < G.vertices.size(); ++v)
    if (G.vertices[v]->n == node->data.root) {
      local_root = v;
      break;
    }
  assert(local_root > -1);
  for (auto H : G.WithoutVertex(local_root))
    reconstruct(H, node->data.root, tree);
}

// Little helper function that returns the treedepth for the given graph.
std::pair<int, std::vector<int>> treedepth(const SubGraph &G) {
  cache = SetTrie();
  int td = treedepth(G, 1, G.vertices.size()).second;
  std::vector<int> tree(G.vertices.size(), -2);
  reconstruct(G, -1, tree);
  // The reconstruction is 0 based, the output is 1 based indexing, fix.
  for (auto &v : tree) v++;
  return {td, std::move(tree)};
}