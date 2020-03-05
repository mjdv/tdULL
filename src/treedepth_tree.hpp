#pragma once

#include <numeric>
#include "graph.hpp"

void root_rank(const SubGraph &G, int v, std::vector<std::set<int>> &L, int d) {
  int a = 0;
  for (auto vi : G.Adj(v)) {
    if (!L[vi].empty()) {
      a = std::max(a, *L[vi].rbegin());
    }
  }
  int c = a + 1;
  int avail = a + 1;
  L[v] = std::set<int>();

  while (G.vertices[v]->rank == -1 && c > -2) {
    c--;

    int i = -1;
    // check in how many lists the max element is equal to c, if only one
    // remember it in i, if more than 1 set i to -2
    for (auto vi : G.Adj(v)) {
      int ti = L[vi].empty() ? 0 : (*L[vi].rbegin());
      if (ti == c) {
        if (i != -1) {
          i = -2;
        } else {
          i = vi;
        }
      }
    }

    if (i == -2 or (c == 0 && d == 1)) {
      G.vertices[v]->rank = avail;

      auto it = L[v].lower_bound(avail);
      L[v].erase(L[v].begin(), it);
      L[v].insert(avail);
    } else if (i > -1) {
      L[v].insert(c);
      L[i].erase(c);
    } else if (i == -1) {
      avail = c;
    }
  }
}

void critical_rank_tree(const SubGraph &G, std::pair<int, int> v,
                        std::vector<std::set<int>> &L) {
  // v.first contains vertex, v.second the previous vertex

  // if you are a leaf node and we're not the first node we're visiting
  if (G.Adj(v.first).size() == 1 && v.first != v.second) {
    G.vertices[v.first]->rank = 1;
    L[v.first] = {1};
  } else {
    int d = 0;
    for (auto vi : G.Adj(v.first)) {
      if (vi != v.second) {
        d++;
        critical_rank_tree(G, std::make_pair(vi, v.first), L);
      }
    }
    root_rank(G, v.first, L, d);
  }
}

std::pair<int, int> treedepth_tree(const SubGraph &G) {
  for (auto v : G.vertices) {
    v->rank = -1;
  }

  // create list for every vertex
  auto L = std::vector<std::set<int>>(G.vertices.size(), std::set<int>());
  critical_rank_tree(G, std::make_pair(0, 0), L);

  // std::cout << "ranks" << std::endl;
  int maxr = 0, root = -1;
  for (auto &v : G.vertices) {
    // std::cout << (v->n+1) << " : " << v->rank << std::endl;
    if (v->rank > maxr) {
      maxr = v->rank;
      root = v->n;
    }
  }
  assert (root > -1);
  return {maxr, root};
}
