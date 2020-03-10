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

  // Add this graph to the cache.
  Node *node;
  bool inserted;
  std::tie(node, inserted) = cache.Insert(G);

  if (!inserted) {
    // This graph was in the cache, retrieve lower/upper bounds.
    lower = node->lower_bound;
    upper = node->upper_bound;
  } else {
    // If this graph wasn't in the cache, store the trivial bounds.
    node->lower_bound = lower;
    node->upper_bound = upper;
    node->root = G.vertices[0]->n;

    // We do a quick check for special cases we can answer exactly and
    // immediately.
    if (N < exactCacheSize) {
      auto [td, root] = exactCache(G.adj);
      return CacheUpdate(node, td, td, G.vertices.at(root)->n);
    } else if (G.IsCompleteGraph()) {
      return CacheUpdate(node, N, N, G.vertices[0]->n);
    } else if (G.IsStarGraph()) {
      // Find node with max_degree.
      for (int v = 0; v < G.vertices.size(); ++v)
        if (G.Adj(v).size() == G.max_degree)
          return CacheUpdate(node, 2, 2, G.vertices[v]->n);
    } else if (G.IsCycleGraph()) {
      // Find the bound, 1 + td of path of length N - 1.
      N--;
      int bnd = 2;
      while (N >>= 1) bnd++;
      return CacheUpdate(node, bnd, bnd, G.vertices[0]->n);
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
          return CacheUpdate(node, bnd, bnd, G.vertices[v]->n);
        }
    } else if (G.IsTreeGraph()) {
      auto [td, root] = treedepth_tree(G);
      return CacheUpdate(node, td, td, root);
    }
  }

  // If the trivial or previously found bounds suffice, we are done.
  if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
    return {lower, upper};
  }

  // If the graph has at least 3 vertices, we never want a leaf (degree 1
  // node) as a root.
  assert(G.vertices.size() > 2);

  std::vector<int> vertices;
  std::vector<int> degree_hist(G.max_degree + 1, 0);
  vertices.reserve(G.vertices.size());
  for (int v = 0; v < G.vertices.size(); ++v) {
    // Only add vertices with deg > 1.
    if (G.Adj(v).size() > 1) vertices.emplace_back(v);
    degree_hist[G.Adj(v).size()]++;
  }
  lower = std::max(lower, (G.M - degree_hist[2]) / (N - degree_hist[2]) + 1);
  assert(vertices.size());

  // Below we have the following two checks.
  // * If G has deg 1 vertices, then we can recursively remove all the deg 1
  // vertices, and call treedepth on the resulting core. This will give a (nice)
  // lower bound.
  // * If G does not have deg 1 vertices, then we can contract all deg 2
  // vertices and try out the trivial lower bound.
  assert(!G.IsTreeGraph());
  for (int k = 2; k <= G.max_degree; k++)
    if (degree_hist[k - 1]) {
      auto cc_core = G.kCore(k);

      // Only try if we do not have an empty core.
      if (!cc_core.empty()) {
        assert(cc_core[0].vertices.size() < N);
        for (const auto &cc : cc_core)
          lower =
              std::max(lower, treedepth(cc, search_lbnd, search_ubnd).first);
      }
      break;
    }

  node->lower_bound = lower;
  if (search_ubnd <= lower || lower == upper) {
    return {lower, upper};
  }

  // If the graph is very dense, it may be that the complement is disconnected.
  // In this case we can reduce the treedepth of G to the treedepths of the
  // subgraphs whose vertex sets are given by the connected components of the
  // complement. (This will be rare, but if it works, very powerful.)
  //
  // In most cases max_degree, N, G.M will rule out that this will work.
  if (G.max_degree * 2 >= N && G.max_degree * (N - G.max_degree) <= G.M) {
    auto components = G.ComplementComponents();
    if (components.size() > 1) {
      // Now td(G) = min_{H : components} |G| - |H| + td(H).
      int lower_components = N;
      for (const auto &H : components) {
        int H_size = H.vertices.size();
        int search_lbnd_H = std::max(1, search_lbnd - N + H_size);
        int search_ubnd_H = std::max(1, search_ubnd - N + H_size);
        auto [lower_H, upper_H] = treedepth(H, search_lbnd_H, search_ubnd_H);

        lower_components = std::min(lower_H + N - H_size, lower_components);
        if (upper_H + N - H_size < upper) {
          upper = upper_H + N - H_size;
          for (int i = 0; i < full_graph_as_sub.vertices.size(); i++) {
            if (G.mask[i] && !H.mask[i]) {
              node->root = i;
              break;
            }
          }
        }
      }
      node->lower_bound = lower = std::max(lower, lower_components);
      node->upper_bound = upper;

      return {lower, upper};
    }
  }

  // Change BetweennessCentrality to DegreeCentrality to go back to the old
  // behaviour of ordering by degree.
  auto centrality = DegreeCentrality(G);

  // Sort the vertices based on the degree.
  std::sort(vertices.begin(), vertices.end(),
            [&](int v1, int v2) { return centrality[v1] > centrality[v2]; });

  if (inserted) {
    // If G is a new graph in the cache, compute its DfsTree-tree from
    // the most promising node once, and then evaluate the treedepth_tree on
    // this tree.
    lower = std::max(lower, treedepth_tree(G.DfsTree(vertices[0])).first);
    node->lower_bound = lower;

    // If the trivial or previously found bounds suffice, we are done.
    if (search_ubnd <= lower || lower == upper) {
      return {lower, upper};
    }
  }

  // Main loop: try every vertex as root.
  // new_lower tries to find a new treedepth lower bound on this subgraph.
  int new_lower = N;
  for (auto v : vertices) {
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
      node->upper_bound = upper;
      node->root = G.vertices[v]->n;
    }

    if (upper <= search_lbnd || lower == upper) {
      // Choosing root v already gives us a treedepth decomposition which is
      // good enough (either a sister branch is at least this long, or it
      // matches a previously proved lower bound for this subgraph) so we
      // can use v as our root.
      return {lower, upper};
    }
  }

  lower = std::max(lower, new_lower);
  node->lower_bound = lower;
  return {lower, upper};
}

// Recursive function to reconstruct the tree that atains the treedepth.
void reconstruct(const SubGraph &G, int root, std::vector<int> &tree, int td) {
  assert(G.vertices.size());
  auto node = cache.Search(G);

  // Ensure that the cache contains the correct node.
  treedepth(G, 1, td);
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
  int td = treedepth(G, 1, G.vertices.size()).second;
  std::vector<int> tree(G.vertices.size(), -2);
  reconstruct(G, -1, tree, td);
  // The reconstruction is 0 based, the output is 1 based indexing, fix.
  for (auto &v : tree) v++;
  return {td, std::move(tree)};
}
