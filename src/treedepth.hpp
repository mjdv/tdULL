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

std::map<std::string, std::map<std::string, int>> treedepth_calls;

// This function returns the treedepth and root of a SubGraph using simple
// heuristics. This only works for special graphs. It returns -1, -1, if
// no such exact heuristic is found.
std::pair<int, int> treedepth_exact(const SubGraph &G) {
  int N = G.vertices.size();
  // We do a quick check for special cases we can answer exactly and
  // immediately.
  if (N < exactCacheSize) {
    treedepth_calls["exact hit"]["small"]++;
    auto [td, root] = exactCache(G.adj);
    return {td, G.vertices[root]->n};
  } else if (G.IsCompleteGraph()) {
    treedepth_calls["exact hit"]["complete"]++;
    return {N, G.vertices[0]->n};
  } else if (G.IsStarGraph()) {
    treedepth_calls["exact hit"]["star"]++;
    // Find node with max_degree.
    for (int v = 0; v < G.vertices.size(); ++v)
      if (G.Adj(v).size() == G.max_degree) return {2, G.vertices[v]->n};
  } else if (G.IsCycleGraph()) {
    treedepth_calls["exact hit"]["cycle"]++;
    // Find the bound, 1 + td of path of length N - 1.
    N--;
    int bnd = 2;
    while (N >>= 1) bnd++;
    return {bnd, G.vertices[0]->n};
  } else if (G.IsPathGraph()) {
    treedepth_calls["exact hit"]["path"]++;
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
        return {bnd, G.vertices[v]->n};
      }
  } else if (G.IsTreeGraph()) {
    treedepth_calls["exact hit"]["tree"]++;
    auto [td, root] = treedepth_tree(G);
    return {td, root};
  }
  return {-1, -1};
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
int max_time_treedepth = 10 * 60;  // A time limit of TEN minuts for now.

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
  if (N > 2 && !G.IsStarGraph()) lower++;
  if (!G.IsCompleteGraph()) upper = N - 1;

  // If the trivial or previously found bounds suffice, we are done.
  if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
    treedepth_calls["calls"]["trivial bounds"]++;
  }

  // Add this graph to the cache.
  Node *node;
  bool inserted;
  std::tie(node, inserted) = cache.Insert(G);

  if (!inserted) {
    // This graph was in the cache, retrieve lower/upper bounds.
    lower = node->lower_bound;
    upper = node->upper_bound;
    treedepth_calls["cache"]["cache hit"]++;
  } else {
    treedepth_calls["cache"]["cache miss"]++;
    // If this graph wasn't in the cache, store the trivial bounds.
    node->lower_bound = lower;
    node->upper_bound = upper;
    node->root = G.vertices[0]->n;

    // Run some checks to see if we can simply find the exact td already.
    auto [td_exact, root_exact] = treedepth_exact(G);
    if (td_exact > -1 && root_exact > -1) {
      treedepth_calls["calls"]["exact hit"]++;
      return CacheUpdate(node, td_exact, td_exact, root_exact);
    }
  }

  // If the trivial or previously found bounds suffice, we are done.
  if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
    treedepth_calls["calls"]["simple bounds"]++;
    return {lower, upper};
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
      node->lower_bound = lower =
          std::max(lower, std::get<0>(treedepth(cc, search_lbnd, search_ubnd)));
      if (search_ubnd <= lower || lower == upper) {
        treedepth_calls["calls"]["kcore recursion"]++;
        return {lower, upper};
      }
    }
  }

  if (inserted) {
    // If the graph has at least 3 vertices, we never want a leaf (degree 1
    // node) as a root.
    assert(G.vertices.size() > 2 && G.max_degree >= 2);
    std::vector<int> vertices;
    vertices.reserve(N);
    for (int v = 0; v < G.vertices.size(); ++v) {
      // Only add vertices with deg > 1.
      if (G.Adj(v).size() > 1) vertices.emplace_back(v);
    }

    // Change BetweennessCentrality to DegreeCentrality to go back to the old
    // behaviour of ordering by degree.
    auto centrality = DegreeCentrality(G);

    // Sort the vertices based on the degree.
    std::sort(vertices.begin(), vertices.end(),
              [&](int v1, int v2) { return centrality[v1] > centrality[v2]; });

    // If G is a new graph in the cache, compute its DfsTree-tree from
    // the most promising node once, and then evaluate the treedepth_tree on
    // this tree.
    for (int v : vertices) {
      node->lower_bound = lower =
          std::max(lower, treedepth_tree(G.DfsTree(v)).first);
    }

    // If the trivial or previously found bounds suffice, we are done.
    if (search_ubnd <= lower || lower == upper) {
      treedepth_calls["calls"]["dfstree"]++;
      return {lower, upper};
    }
  }

  // Main loop: try every separator as a set of roots.
  // new_lower tries to find a new treedepth lower bound on this subgraph.
  int new_lower = N;

  for (auto separator : G.AllMinimalSeparators()) {
    int sep_size = separator.size();

    int search_ubnd_sep =
        std::max(1, std::min(search_ubnd - sep_size, upper - sep_size));
    int search_lbnd_sep = std::max(search_lbnd - sep_size, 1);

    int upper_sep = 0;
    int lower_sep = lower - sep_size;

    bool early_break = false;
    for (auto H : G.WithoutVertices(separator)) {
      auto [lower_H, upper_H] = treedepth(H, search_lbnd_sep, search_ubnd_sep);

      upper_sep = std::max(upper_sep, upper_H);
      lower_sep = std::max(lower_sep, lower_H);

      search_lbnd_sep = std::max(search_lbnd_sep, lower_H);

      if (lower_H > search_ubnd_sep) {
        // This component already shows that there's no reason to
        // continue trying with vertex v.
        early_break = true;
        break;
      }
    }

    new_lower = std::min(new_lower, lower_sep + sep_size);

    // If we find a new upper bound, update the cache accordingly :-).
    if (!early_break && upper_sep + sep_size < upper) {
      node->upper_bound = upper = upper_sep + sep_size;
      node->root = G.vertices[separator[0]]->n;

      // Iteratively remove the seperator from G and update bounds.
      SubGraph H = G;
      for (int i = 1; i < separator.size(); i++) {
        // Get the subgraph after removing seperator[i-1].
        auto cc = H.WithoutVertex(H.LocalIndex(G.vertices[separator[i - 1]]));
        if (cc.size() > 1) break;
        H = cc[0];
        auto [node_H, inserted_H] = cache.Insert(H);

        // Now if H was new to the cache, or we have better bounds, lets update!
        if (inserted_H || (upper - i < node_H->upper_bound)) {
          // Check if we can find tight bounds for H directly.
          auto [td_exact_H, root_exact_H] = treedepth_exact(H);
          if (td_exact_H > -1 && root_exact_H > -1)
            CacheUpdate(node_H, td_exact_H, td_exact_H, root_exact_H);
          else {
            // No exact bounds, so simply store our lower/upper bounds.
            node_H->upper_bound = upper - i;
            node_H->lower_bound = std::max(lower - i, node_H->lower_bound);
            node_H->root = G.vertices[separator[i]]->n;
          }
        }
      }
    }

    if (upper <= search_lbnd || lower == upper) {
      // Choosing root v already gives us a treedepth decomposition which is
      // good enough (either a sister branch is at least this long, or it
      // matches a previously proved lower bound for this subgraph) so we
      // can use v as our root.
      treedepth_calls["calls"]["separators recursion"]++;
      return {lower, upper};
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
  return {lower, upper};
}

// Recursive function to reconstruct the tree that atains the treedepth.
void reconstruct(const SubGraph &G, int root, std::vector<int> &tree, int td) {
  assert(G.vertices.size());
  auto node = cache.Search(G);

  // Ensure that the cache contains the correct node.
  treedepth(G, td, G.vertices.size());
  node = cache.Search(G);
  assert(node);
  assert(node->root > -1);
  tree.at(node->root) = root;

  // Root is the global coordinate, find its local coordinate.
  int local_root = -1;
  for (int v = 0; v < G.vertices.size(); ++v)
    if (G.vertices[v]->n == node->root) {
      local_root = v;
      break;
    }
  assert(local_root > -1);
  for (auto H : G.WithoutVertex(local_root))
    reconstruct(H, node->root, tree, td - 1);
}

// Little helper function that returns the treedepth for the given graph.
std::pair<int, std::vector<int>> treedepth(const SubGraph &G) {
  cache = SetTrie();
  treedepth_calls.clear();

  time(&time_start_treedepth);
  int td = treedepth(G, 1, G.vertices.size()).second;
  std::vector<int> tree(G.vertices.size(), -2);
  reconstruct(G, -1, tree, td);
  // The reconstruction is 0 based, the output is 1 based indexing, fix.
  for (auto &v : tree) v++;

  treedepth_calls["cache"]["size"] = cache.size();
  std::cout << "Graph recursion data: " << std::endl;
  for (const auto &[key, map] : treedepth_calls) {
    std::cout << "\t" << key << ": ";
    for (const auto &[call, value] : map) {
      std::cout << call << " - " << value << "; ";
    }
    std::cout << std::endl;
  }
  return {td, std::move(tree)};
}
