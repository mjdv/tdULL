#include <array>
#include <bitset>
#include <iostream>
#include <memory>
#include <vector>

#include "treedepth.hpp"

using ull = unsigned long long;
constexpr ull N = 8;
constexpr ull M = N * (N - 1) / 2;
std::string filename = "exact_cache_" + std::to_string(N) + ".bin";

std::array<std::array<int16_t, N>, N> mapping;

bool isConnected(std::bitset<M> &edges) {
  std::vector<int> stack;
  std::vector<bool> visited;
  stack.reserve(N);
  visited.resize(N, false);

  stack.push_back(0);
  visited[0] = true;

  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();

    for (int w = 0; w < N; ++w) {
      if (v == w) continue;
      int x = std::min(v, w);
      int y = std::max(v, w);
      ull edge = mapping[x][y];
      assert(edge >= 0 && edge < M);
      if (edges[edge] && !visited[w]) {
        visited[w] = true;
        stack.push_back(w);
      }
    }
  }

  for (int v = 0; v < N; ++v)
    if (!visited[v]) return false;

  return true;
}

// TODO Ray check of dit nog klopt
int main() {
  std::vector<int> vertices;
  for (int v = 0; v < N; v++) vertices.emplace_back(v);

  for (int v = 0; v < N; ++v)
    for (int w = 0; w < N; ++w) mapping[v][w] = -1;

  int num = 0;
  for (int v = 0; v < N; ++v)
    for (int w = v + 1; w < N; ++w) mapping[v][w] = num++;

  ull total = 0;
  ull totalcc = 0;
  time_t start, now;
  time(&start);
  std::ofstream output(filename, std::ios::binary);
  for (ull i = 0; i < (1 << M); ++i) {
    std::bitset<M> edges(i);
    total++;
    if (isConnected(edges)) {
      totalcc++;
      Graph sub;
      sub.global = vertices;
      sub.M = edges.count();
      sub.adj.resize(N);
      for (int v = 0; v < N; ++v) {
        for (int w = v + 1; w < N; ++w) {
          ull edge = mapping[v][w];
          if (edges[edge]) {
            sub.adj[v].push_back(w);
            sub.adj[w].push_back(v);
          }
        }
        sub.max_degree = std::max(sub.max_degree, sub.adj[v].size());
        sub.min_degree = std::min(sub.min_degree, sub.adj[v].size());
      }
      auto [td, _] = treedepth(sub);
      int root = cache.Search(sub)->root;
      assert(root > -1 && td >= 1 && td <= N);
      output.put(uint8_t((td << 4) | root));
    } else {
      output.put(0);
    }
    if (i % 10000 == 0) {
      time(&now);
      std::cout << "Total connected graphs " << totalcc << " / " << i
                << ". 10k graphs/s:  "
                << difftime(now, start) / double(i / 10000) << "." << std::endl;
    }
  }

  std::cout << "Total connected graphs " << totalcc << " / " << total
            << std::endl;
  return 0;
}
