#pragma once
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>

#include "centrality.hpp"
#include "exact_cache.hpp"
#include "graph.hpp"
#include "set_trie.hpp"
#include "treedepth_tree.hpp"

// Trivial treedepth implementation, useful for simple sanity checks.
int treedepth_trivial(const SubGraph &G) {
  if (G.vertices.size() == 1) return 1;
  int td = G.vertices.size();
  for (int v = 0; v < G.vertices.size(); ++v) {
    int td_for_this_root = 0;
    for (auto cc : G.WithoutVertex(v)) {
      td_for_this_root = std::max(td_for_this_root, treedepth_trivial(cc) + 1);
    }
    td = std::min(td, td_for_this_root);
  }
  return td;
}

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

// Little helper function to update information in the cache.
std::pair<int, int> CacheUpdate(Node *node, int lower_bound, int upper_bound,
                                int root) {
  node->lower_bound = lower_bound;
  node->upper_bound = upper_bound;
  node->root = root;
  return std::pair{lower_bound, upper_bound};
};

// Global variable keeping track of the time we've spent so far, and the limit.
time_t time_start_treedepth;
int max_time_treedepth = 10 * 60;  // Lets put the time limit it at ten minutes for now.

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
std::tuple<int, int, int> treedepth(const SubGraph &G, int search_lbnd,
                                    int search_ubnd) {
  int N = G.vertices.size();
  if (N == 1) return {1, 1, G.vertices[0]->n};
  if (search_lbnd >= search_ubnd) return {G.M / N + 1, N, G.vertices[0]->n};

  // We do a quick check for special cases we can answer exactly and
  // immediately.
  if (N < exactCacheSize) {
    auto [td, root] = exactCache(G.adj);
    return {td, td, G.vertices.at(root)->n};
  } else if (G.IsCompleteGraph()) {
    return {N, N, G.vertices[0]->n};
  } else if (G.IsStarGraph()) {
    // Find node with max_degree.
    for (int v = 0; v < G.vertices.size(); ++v)
      if (G.Adj(v).size() == G.max_degree) return {2, 2, G.vertices[v]->n};
  } else if (G.IsCycleGraph()) {
    // Find the bound, 1 + td of path of length N - 1.
    N--;
    int bnd = 2;
    while (N >>= 1) bnd++;
    return {bnd, bnd, G.vertices[0]->n};
  } else if (G.IsPathGraph()) {
    // Find the bound, this is the ceil(log_2(N)).
    int bnd = 1;
    while (N >>= 1) bnd++;

    // Find a leaf and then find the middle node.
    for (int v = 0; v < G.vertices.size(); ++v)
      if (G.Adj(v).size() == 1) {
        int prev = v;
        v = G.adj[v][0];

        // Find the middle node.
        for (int i = 1; i < G.vertices.size() / 2; i++) {
          int tmp = v;
          v = (prev ^ G.adj[v][0] ^ G.adj[v][1]);
          prev = tmp;
        }
        return {bnd, bnd, G.vertices[v]->n};
      }
  } else if (G.IsTreeGraph()) {
    // TODO: this one is semi-expensive, but probably doesn't occur often.
    auto [td, root] = treedepth_tree(G);
    return {td, td, root};
  }

  // None of the exact tests worked, so lets try the trivial bounds.
  int lower = G.M / N + 1;
  int upper = N - 1;
  int root = G.vertices[0]->n;

  // If the trivial or previously found bounds suffice, we are done.
  if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
    return {lower, upper, root};
  }

  // If not, lets check if it already exists in the cache.
  Node *node = cache.Search(G);
  if (node) {
    // This graph was in the cache, retrieve lower/upper bounds.
    lower = node->lower_bound;
    upper = node->upper_bound;
    root = node->root;

    if (search_ubnd <= lower || search_lbnd >= upper || lower == upper)
      return {lower, upper, root};
  }

  // Below we calculate the smallest k-core that G can contain. If this is non-
  // empty, we recursively calculate the treedepth on this core first. This
  // should give a nice lower bound pretty rapidly.
  auto cc_core = G.kCore(G.min_degree + 1);
  if (!cc_core.empty()) {
    assert(cc_core[0].vertices.size() < N);

    // Sort the components on density.
    std::sort(cc_core.begin(), cc_core.end(), [](auto &c1, auto &c2) {
      return c1.M / c1.vertices.size() > c2.M / c2.vertices.size();
    });
    for (const auto &cc : cc_core) {
      lower =
          std::max(lower, std::get<0>(treedepth(cc, search_lbnd, search_ubnd)));
      if (search_ubnd <= lower || lower == upper) return {lower, upper, root};
    }
  }

  // If the graph has at least 3 vertices, we never want a leaf (degree 1
  // node) as a root.
  assert(G.vertices.size() > 2 && G.max_degree >= 2);
  std::vector<int> vertices;
  vertices.reserve(N);
  int N_deg_2 = 0;
  for (int v = 0; v < G.vertices.size(); ++v) {
    // Only add vertices with deg > 1.
    if (G.Adj(v).size() > 1) vertices.emplace_back(v);

    // Count number of deg 2 vertices.
    if (G.Adj(v).size() == 2) N_deg_2++;
  }
  lower = std::max(lower, (G.M - N_deg_2) / (N - N_deg_2) + 1);
  if (search_ubnd <= lower || lower == upper) return {lower, upper, root};

  // Change DegreeCentrality to other centralities here, say BetweennessCentrality.
  auto centrality = DegreeCentrality(G);

  // Sort the vertices based on the degree.
  std::sort(vertices.begin(), vertices.end(),
            [&](int v1, int v2) { return centrality[v1] > centrality[v2]; });

  // If G doesn't exist in the cache, lets add it now, since we will start doing
  // some real work.
  if (node == nullptr) {
    node = cache.Insert(G).first;
    node->lower_bound = lower;
    node->upper_bound = upper;
    node->root = G.vertices[0]->n;

    // If G is a new graph in the cache, compute its DfsTree-tree from
    // the most promising node once, and then evaluate the treedepth_tree on
    // this tree.
    node->lower_bound = lower =
        std::max(lower, treedepth_tree(G.DfsTree(vertices[0])).first);
    if (search_ubnd <= lower || lower == upper) return {lower, upper, root};
  }

  // Main loop: try every vertex as root.
  // new_lower tries to find a new treedepth lower bound on this subgraph.
  int new_lower = N;
  for (auto v : vertices) {
    int search_ubnd_v = std::min(search_ubnd - 1, upper - 1);
    int search_lbnd_v = std::max(search_lbnd - 1, 1);

    int upper_v = 0;
    int lower_v = lower - 1;

    // Sort the components of G \ {v} on density.
    auto cc = G.WithoutVertex(v);
    std::sort(cc.begin(), cc.end(), [](const SubGraph &c1, const SubGraph &c2) {
      return c1.M / c1.vertices.size() > c2.M / c2.vertices.size();
    });

    // Recursively find the treedepths for each of the component.
    for (const auto &H : cc) {
      auto [lower_H, upper_H, _] =
          treedepth(H, search_lbnd_v, search_ubnd_v);

      upper_v = std::max(upper_v, upper_H);
      lower_v = std::max(lower_v, lower_H);

      search_lbnd_v = std::max(search_lbnd_v, lower_H);
    }

    new_lower = std::min(new_lower, lower_v + 1);

    // If we found a new upper bound, update the root accordingly.
    if (upper_v + 1 < upper) {
      node->upper_bound = upper = upper_v + 1;
      node->root = root = G.vertices[v]->n;
    }

    if (upper <= search_lbnd || lower == upper) {
      // Choosing root v already gives us a treedepth decomposition which is
      // good enough (either a sister branch is at least this long, or it
      // matches a previously proved lower bound for this subgraph) so we
      // can use v as our root.
      return {lower, upper, root};
    }

    // Check whether we are still in the time limits.
    time_t now;
    time(&now);
    if (difftime(now, time_start_treedepth) > max_time_treedepth)
      throw std::runtime_error(
          "Ran out of time, spent " +
          std::to_string(difftime(now, time_start_treedepth)) + " seconds.");
  }

  node->lower_bound = lower = std::max(lower, new_lower);
  return {lower, upper, root};
}

// Recursive function to reconstruct the tree that atains the treedepth.
void reconstruct(const SubGraph &G, int root, std::vector<int> &tree, int td) {
  assert(G.vertices.size());

  // Ensure that the cache contains the correct node.
  int new_root = std::get<2>(treedepth(G, 1, td));
  assert(new_root > -1);
  tree.at(new_root) = root;

  // Root is the global coordinate, find its local coordinate.
  int local_root = -1;
  for (int v = 0; v < G.vertices.size(); ++v)
    if (G.vertices[v]->n == new_root) {
      local_root = v;
      break;
    }
  assert(local_root > -1);
  for (auto H : G.WithoutVertex(local_root))
    reconstruct(H, new_root, tree, td - 1);
}

// Little helper function that returns the treedepth for the given graph.
std::pair<int, std::vector<int>> treedepth(const SubGraph &G) {
  cache = SetTrie();
  time(&time_start_treedepth);
  int td = std::get<1>(treedepth(G, 1, G.vertices.size()));
  std::vector<int> tree(G.vertices.size(), -2);
  reconstruct(G, -1, tree, td);
  // The reconstruction is 0 based, the output is 1 based indexing, fix.
  for (auto &v : tree) v++;
  return {td, std::move(tree)};
}
