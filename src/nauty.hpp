#pragma once
#include "graph.hpp"

extern "C" {
#include "../third_party/nauty/nausparse.h"
}

struct group_struct;
struct Nauty {
  Nauty(const Graph &G, bool canonical_labeling = false);
  ~Nauty();
  const std::vector<std::vector<int>> &Automorphisms();

  // Stats.
  size_t num_automorphisms = 0;
  size_t num_generators = 0;
  size_t num_orbits = 0;

  // Calculated values.
  std::vector<int> orbit_representatives;
  std::vector<int> orbits;
  std::vector<std::vector<int>> automorphisms;
  long canon_hash = 0;
  sparsegraph cg;

 protected:
  const Graph &G;
  group_struct *group = nullptr;
};

struct NautyHashData {
  sparsegraph cg;  // Graph with canonical labeling.
  int td = -1;
  int root = -1;

  NautyHashData(sparsegraph *nauty_sg) {
    SG_INIT(cg);
    copy_sg(nauty_sg, &cg);
  }
};

// Simple graph cache.
class NautyHashCache {
 public:
  // Tries to find the given graph, returns pointer to data if it succeeded.
  NautyHashData *Search(Nauty &nauty) {
    assert(nauty.canon_hash != 0);
    auto it = cache_.find(nauty.canon_hash);

    // If data doesn't yet exist in cache, we can stop looking.
    if (it == cache_.end()) return nullptr;

    NautyHashData &data = it->second;

    // Check if the stored graphs coincide.
    assert(aresame_sg(&nauty.cg, &data.cg));

    return &data;
  }

  // Inserts data into the cache if it did not yet exists,
  // returns pointer if it succesfully inserted it.
  NautyHashData *Insert(Nauty &nauty) {
    assert(nauty.canon_hash != 0);
    auto it = cache_.find(nauty.canon_hash);

    // If data doesn't yet exist in cache, simply add it.
    if (it == cache_.end())
      return &cache_.emplace(nauty.canon_hash, NautyHashData(&nauty.cg))
                  .first->second;
    return nullptr;
  }

 protected:
  std::map<long, NautyHashData> cache_;
};
