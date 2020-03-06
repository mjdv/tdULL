#include "graph_hash.hpp"

#include <algorithm>
#include <iostream>
#include <random>

std::vector<std::vector<int>> RandomPermute(
    const std::vector<std::vector<int>>& G) {
  std::vector<std::vector<int>> result(G.size());

  std::vector<int> mapping(G.size());
  std::iota(mapping.begin(), mapping.end(), 0);

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(mapping.begin(), mapping.end(), g);

  for (int v_old = 0; v_old < G.size(); ++v_old) {
    int v_new = mapping[v_old];
    std::vector<int> nghbrs;
    for (int nghb : G[v_old]) nghbrs.push_back(mapping[nghb]);
    result[v_new] = nghbrs;
  }
  return result;
}

std::vector<std::vector<int>> RandomGraph(int N, double p) {
  std::vector<std::vector<int>> G(N);
  for (int v = 0; v < N; ++v)
    for (int w = v + 1; w < N; ++w) {
      double r = rand() / double(RAND_MAX);
      if (r <= p) {
        G[v].push_back(w);
        G[w].push_back(v);
      }
    }
  return G;
}

bool isConnected(const std::vector<std::vector<int>>& adj) {
  int N = adj.size();
  std::vector<int> stack;
  std::vector<bool> visited;
  stack.reserve(N);
  visited.resize(N, false);

  stack.push_back(0);
  visited[0] = true;

  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();

    for (auto w : adj[v]) {
      if (!visited[w]) {
        visited[w] = true;
        stack.push_back(w);
      }
    }
  }

  for (int v = 0; v < N; ++v)
    if (!visited[v]) return false;

  return true;
}

bool isCycle(const std::vector<std::vector<int>>& adj) {
  // Assume adj is connected.
  for (auto bla : adj)
    if (bla.size() != 2) return false;
  return true;
}

bool isPath(const std::vector<std::vector<int>>& adj) {
  int M = 0;
  int N = adj.size();
  // Assume adj is connected.
  for (auto bla : adj) {
    M += bla.size();
    if (bla.size() > 2) return false;
  }
  M /= 2;
  return N - 1 == M;
}

std::ostream& operator<<(std::ostream& os,
                         const std::vector<std::vector<int>>& adj) {
  for (int i = 0; i < adj.size(); ++i) {
    os << i << " : {";
    for (int j = 0; j < adj[i].size(); ++j) {
      if (j) os << ", ";
      os << adj[i][j];
    }
    os << "}" << std::endl;
  }
  return os;
}

int main() {
  std::vector<std::vector<int>> triangle{{1, 2}, {0}, {0}};
  for (int i = 0; i < 10; ++i) {
    auto permuted = RandomPermute(triangle);
    assert(GraphHashIsomorphism(triangle, permuted));
  }

  srand(0);
  int iso_misses_graphs = 0;
  int iso_misses_total = 0;
  int hash_misses = 0;
  std::vector<std::vector<int>> G;
  for (int j = 0; j < 1000; ++j) {
    int N = 10;

    /// Find a random connected graph that is not a circle.
    while (true) {
      G = RandomGraph(N, 0.5);
      if (isConnected(G) && !isCycle(G) && !isPath(G)) break;
    }
    int M = 0;
    for (auto adj : G) M += adj.size();
    assert(M % 2 == 0);
    M /= 2;

    assert(GraphHashIsomorphism(G, G));

    bool missed = false;
    auto G_hashes = GraphHash(G);
    for (int i = 0; i < 50; ++i) {
      auto permuted = RandomPermute(G);
      auto permuted_hashes = GraphHash(permuted);
      if (G_hashes.first != permuted_hashes.first) hash_misses++;

      // Try multiple times to find a correct isomphism mapping.
      bool found_isomophism = false;
      for (int m = 0; m < 10; ++m) {
        auto [succes, mapping] = GraphHashIsomphismMapping(
            G, permuted, G_hashes.second, permuted_hashes.second,
            /* randomize */ true);

        // Assert that we can find _a_ mapping, might not be correct.
        assert(succes);

        if (GraphIsomorphism(G, permuted, mapping)) {
          found_isomophism = true;
          break;
        }
      }

      if (!found_isomophism) {
        iso_misses_total++;
        missed = true;
      }
    }
    if (missed) iso_misses_graphs++;
  }
  std::cout << "Total of " << hash_misses
            << " hash misses for isomoprhic graphs." << std::endl;
  std::cout << "Total of " << iso_misses_total << " isomoprhism misses for "
            << iso_misses_graphs << " different graphs" << std::endl;
}
