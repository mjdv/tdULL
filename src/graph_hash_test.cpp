#include "graph_hash.hpp"

#include <algorithm>
#include <iostream>
#include <random>

std::vector<std::vector<int>> RandomPermute(
    const std::vector<std::vector<int>> &G) {
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

int main() {
  std::vector<std::vector<int>> triangle{{1, 2}, {0}, {0}};
  for (int i = 0; i < 10; ++i) {
    auto permuted = RandomPermute(triangle);
    assert(GraphHashIsomorphism(triangle, permuted));
  }

  srand(0);
  int misses_total = 0;
  int misses_per_g = 0;
  for (int j = 0; j < 1000; ++j) {
    int N = 4;
    auto G = RandomGraph(N, 0.5);
    int M = 0;
    for (auto adj : G) M += adj.size();
    assert(M % 2 == 0);
    M /= 2;

    assert(GraphHashIsomorphism(G, G));

    bool missed = false;
    for (int i = 0; i < 50; ++i) {
      auto permuted = RandomPermute(G);
      if (!GraphHashIsomorphism(G, permuted)) {
        misses_total++;
        missed = true;
      }
    }
    if (missed) misses_per_g++;
  }
  std::cout << misses_total << " " << misses_per_g << std::endl;
}
