#include "exact_cache.hpp"

#include "exact_caches/incbin.h"

#if EXACT_CACHE_SIZE >= 2
INCBIN(ExactCache2, "exact_caches/exact_cache_2.bin");
#endif
#if EXACT_CACHE_SIZE >= 3
INCBIN(ExactCache3, "exact_caches/exact_cache_3.bin");
#endif
#if EXACT_CACHE_SIZE >= 4
INCBIN(ExactCache4, "exact_caches/exact_cache_4.bin");
#endif
#if EXACT_CACHE_SIZE >= 5
INCBIN(ExactCache5, "exact_caches/exact_cache_5.bin");
#endif
#if EXACT_CACHE_SIZE >= 6
INCBIN(ExactCache6, "exact_caches/exact_cache_6.bin");
#endif
#if EXACT_CACHE_SIZE >= 7
INCBIN(ExactCache7, "exact_caches/exact_cache_7.bin");
#endif
#if EXACT_CACHE_SIZE >= 8
INCBIN(ExactCache8, "exact_caches/exact_cache_8.bin");
#endif

const std::vector<std::vector<int>> &exactCacheMapping(int N) {
  static std::array<std::vector<std::vector<int>>, exactCacheSize> mappings;
  assert(N < EXACT_CACHE_SIZE);
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
      return {1, 0};
#if EXACT_CACHE_SIZE >= 2
    case 2:
      assert(edges < gExactCache2Size);
      result = gExactCache2Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 3
    case 3:
      assert(edges < gExactCache3Size);
      result = gExactCache3Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 4
    case 4:
      assert(edges < gExactCache4Size);
      result = gExactCache4Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 5
    case 5:
      assert(edges < gExactCache5Size);
      result = gExactCache5Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 6
    case 6:
      assert(edges < gExactCache6Size);
      result = gExactCache6Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 7
    case 7:
      assert(edges < gExactCache7Size);
      result = gExactCache7Data[edges];
      break;
#endif
#if EXACT_CACHE_SIZE >= 8
    case 8:
      assert(edges < gExactCache8Size);
      result = gExactCache8Data[edges];
      break;
#endif
  }
  int td = result >> 4;
  int root = result & 15;
  return {td, root};
}
