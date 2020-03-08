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

int main() {
  full_graph.N = N;
  for (int v = 0; v < N; v++) full_graph.vertices.emplace_back(v);
  std::vector<Vertex *> vertices;
  for (int v = 0; v < N; v++) vertices.emplace_back(&full_graph.vertices[v]);

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
      SubGraph sub;
      sub.vertices = vertices;
      sub.M = edges.count();
      sub.mask = std::vector<bool>(N, true);
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
      }
      auto [td, _] = treedepth(sub);
      int root = cache.Search(sub)->data.root;
      // std::cerr << ((td << 4) | root) << ",";
      output.put(uint8_t((td << 4) | root));
    } else {
      // std::cerr << 0 << ",";
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
