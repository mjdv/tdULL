#include "exact_cache.hpp"

uint8_t exact_cache_2[] = {
#include "exact_cache_2.ipp"
};
uint8_t exact_cache_3[] = {
#include "exact_cache_3.ipp"
};
uint8_t exact_cache_4[] = {
#include "exact_cache_4.ipp"
};
uint8_t exact_cache_5[] = {
#include "exact_cache_5.ipp"
};
uint8_t exact_cache_6[] = {
#include "exact_cache_6.ipp"
};
uint8_t exact_cache_7[] = {
#include "exact_cache_7.ipp"
};

const std::vector<std::vector<int>> &exactCacheMapping(int N) {
  static std::array<std::vector<std::vector<int>>, exactCacheSize> mappings;
  assert(N < exactCacheSize);
  if (mappings[N].empty()) {
    mappings[N] = std::vector<std::vector<int>>(N, std::vector<int>(N, -1));
    int num = 0;
    for (int v = 0; v < N; ++v)
      for (int w = v + 1; w < N; ++w) mappings[N][v][w] = num++;
  }
  return mappings[N];
}

std::pair<int, int> exactCache(const std::vector<std::vector<int>> &adj) {
  int N = adj.size();
  assert(N >= 0 && N < exactCacheSize);
  int32_t edges = 0;
  auto &mapping = exactCacheMapping(N);
  for (int v = 0; v < N; ++v) {
    for (int w : adj[v]) {
      if (v < w) {
        int edge = mapping[v][w];
        edges |= (1 << edge);
      }
    }
  }

  uint8_t result = 0;
  switch (N) {
    case 0:
      return {0, 0};
    case 1:
      return {0, 0};
    case 2:
      assert(edges < sizeof(exact_cache_2));
      result = exact_cache_2[edges];
      break;
    case 3:
      assert(edges < sizeof(exact_cache_3));
      result = exact_cache_3[edges];
      break;
    case 4:
      assert(edges < sizeof(exact_cache_4));
      result = exact_cache_4[edges];
      break;
    case 5:
      assert(edges < sizeof(exact_cache_5));
      result = exact_cache_5[edges];
      break;
    case 6:
      assert(edges < sizeof(exact_cache_6));
      result = exact_cache_6[edges];
      break;
    case 7:
      assert(edges < sizeof(exact_cache_7));
      result = exact_cache_7[edges];
      break;
  }
  int td = result >> 4;
  int root = result & 15;
  return {td, root};
}
