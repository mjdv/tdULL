#pragma once

#include <numeric>
#include "graph.hpp"

/*
 *  Implementation of the algorithm described in:
 *  Optimal node ranking of trees
 *  Ananth V. Iyer, H. Donald Ratliff, G. Vijayan
 *  https://doi.org/10.1016/0020-0190(88)90194-9
 *
 *  Runtime should be O(n log(n))
 *
 *  Possible changes/improvements:
 *  - Use a 'local' list to remember ranks, instead of having a field inside a vertex.
 *  - Do some bookkeeping to remember the node with highest rank. 
 *  - Are the right datastructures used, especially for L? 
 *  - Generate the critical rank lists L on the fly?
 *
 *  - Use the linear time algorithm as described in
 *    Optimal node ranking of trees in linear time
 *    Alejandro A. Schaffer
 *    https://doi.org/10.1016/0020-0190(89)90161-0
 *
 *
 */


void root_rank(const Graph &G, int v, int v_prev, std::vector<std::set<int>> &L, int d,
    std::vector<int> &rank) {
  int a = 0;
  // a will contain the maximum over all elements of the L of neighbors of v
  for (auto vi : G.Adj(v)) {
    if (vi != v_prev and !L[vi].empty()) {
      a = std::max(a, *L[vi].rbegin());
    }
  }
  int c = a + 1;
  int avail = a + 1;
  L[v] = std::set<int>();

  while (rank[v] == -1) {
    c--;

    int i = -1;
    // check in how many lists the max element is equal to c, if only one
    // remember it in i, if more than 1 set i to -2
    for (auto vi : G.Adj(v)) {
      if(vi == v_prev) continue;
      int ti = L[vi].empty() ? 0 : (*L[vi].rbegin());
      if (ti == c) {
        if (i != -1) {
          i = -2;
        } else {
          i = vi;
        }
      }
    }

    if (i == -2 or (c == 0 and d == 1)) {
      // if there are two subtrees with critical rank c, then v needs to get a higher rank
      rank[v] = avail;

      auto it = L[v].lower_bound(avail);
      L[v].erase(L[v].begin(), it);
      L[v].insert(avail);
    } else if (i > -1) {
      // if there is only one subtree with the current value of c, we can go lower
      L[v].insert(c);
      L[i].erase(c);
    } else if (i == -1) {
      // No subtree has a rank of c, we can lower the possible rank we're going to assign
      avail = c;
    }
  }
}

void critical_rank_tree(const Graph &G, int v, int v_prev, 
    std::vector<std::set<int>> &L, std::vector<int> &rank) {

  // if you are a leaf node and we're not the first node we're visiting
  if (G.Adj(v).size() == 1 and v != v_prev) {
    rank[v] = 1;
    L[v] = {1};
  } else {
    int d = 0;
    for (auto vi : G.Adj(v)) {
      // for every subtree starting at a neighbor of v, compute the ranking + critical rank list
      if (vi != v_prev) {
        d++;
        critical_rank_tree(G, vi, v, L, rank);
      }
    }
    root_rank(G, v, v_prev, L, d, rank);
  }
}

std::pair<int, int> treedepth_tree(const Graph &G) {
  // Create vector of ranks, that we will pass around.
  int N = G.N;
  std::vector<int> rank(N, -1);

  // Create list for every vertex
  auto L = std::vector<std::set<int>>(N, std::set<int>());

  // Start from vertex 0, can start from any vertex
  critical_rank_tree(G, 0, 0, L, rank);

  // Find the root with the max rank
  int maxr = 0, root = -1;
  for (int v = 0; v < N; v++) {
    if (rank[v] > maxr) {
      maxr = rank[v];
      root = G.global[v];
    }
  }
  assert(root > -1);
  return {maxr, root};
}
