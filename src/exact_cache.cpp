#include "exact_cache.hpp"

#include "exact_caches/incbin.h"
INCBIN(ExactCache2, "exact_caches/exact_cache_2.bin");
INCBIN(ExactCache3, "exact_caches/exact_cache_3.bin");
INCBIN(ExactCache4, "exact_caches/exact_cache_4.bin");
INCBIN(ExactCache5, "exact_caches/exact_cache_5.bin");
INCBIN(ExactCache6, "exact_caches/exact_cache_6.bin");
INCBIN(ExactCache7, "exact_caches/exact_cache_7.bin");
INCBIN(ExactCache8, "exact_caches/exact_cache_8.bin");

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
      assert(edges < gExactCache2Size);
      result = gExactCache2Data[edges];
      break;
    case 3:
      assert(edges < gExactCache3Size);
      result = gExactCache3Data[edges];
      break;
    case 4:
      assert(edges < gExactCache4Size);
      result = gExactCache4Data[edges];
      break;
    case 5:
      assert(edges < gExactCache5Size);
      result = gExactCache5Data[edges];
      break;
    case 6:
      assert(edges < gExactCache6Size);
      result = gExactCache6Data[edges];
      break;
    case 7:
      assert(edges < gExactCache7Size);
      result = gExactCache7Data[edges];
      break;
    case 8:
      assert(edges < gExactCache8Size);
      result = gExactCache8Data[edges];
      break;
  }
  int td = result >> 4;
  int root = result & 15;
  return {td, root};
}
