#pragma once
#include <cstdint>
#include <map>
#include <vector>

// Calculates the graph hash for graph G given by its adj lists.
// Evaluates a maximum of N steps of hash updating.
std::pair<uint32_t, std::vector<uint32_t>> GraphHash(
    const std::vector<std::vector<int>> &G);

// Heuristically calculates a isomoprhism between two graphs on basis of the
// given vertex hashes, if randomize is set to true, this will return
// different mapping in case there are multiple nodes with the same
// hash.
//
// Returns whether it found _a_ mapping, and the given mapping in case it did.
std::pair<bool, std::vector<int>> GraphHashIsomphismMapping(
    const std::vector<std::vector<int>> &G1,
    const std::vector<std::vector<int>> &G2,
    const std::vector<uint32_t> &vertex_h1,
    const std::vector<uint32_t> &vertex_h2, bool randomize = false);

// Calculates whether two graphs are an isomoprhism on basis of the given
// mapping.
bool GraphIsomorphism(const std::vector<std::vector<int>> &G1,
                      const std::vector<std::vector<int>> &G2,
                      const std::vector<int> &mapping);

// Little helper function for convenience.
bool GraphHashIsomorphism(const std::vector<std::vector<int>> &G1,
                          const std::vector<std::vector<int>> &G2);

struct GraphHashData {
  int td = -1;
  int root = -1;
};

// Simple graph hash cache.
class GraphHashCache {
 public:
  // Tries to find the given graph, returns pointer to data and
  // isomoprhism mapping if it succeeded.
  std::tuple<GraphHashData *, std::vector<int>> Search(
      const std::vector<std::vector<int>> &G, int random_tries = 5) {
    auto [G_hash, vertex_hashes] = GraphHash(G);
    auto it = cache_.find(G_hash);

    // If data doesn't yet exist in cache, we can stop looking.
    if (it == cache_.end()) return {nullptr, {}};

    // Else try to match the cached data to the hash.
    auto &[G_match, G_match_hashes, data] = it->second;

    for (int m = 0; m < random_tries; ++m) {
      auto [succes, mapping] =
          GraphHashIsomphismMapping(G, G_match, vertex_hashes, G_match_hashes,
                                    /* randomize */ true);

      if (!succes) continue;

      // We found an isomoprhism with the cached data!
      if (GraphIsomorphism(G, G_match, mapping)) return {&data, mapping};
    }

    // The graph hash exists in the cache, but we couldn't find an isomoprhism.
    return {nullptr, {}};
  }

  // Inserts data into the cache if it did not yet exists,
  // returns pointer if it succesfully inserted it.
  GraphHashData *Insert(const std::vector<std::vector<int>> &G) {
    auto [G_hash, vertex_hashes] = GraphHash(G);
    auto it = cache_.find(G_hash);

    // If data doesn't yet exist in cache, simply add it.
    if (it == cache_.end())
      return &std::get<2>(
          cache_.emplace(G_hash, std::tuple{G, vertex_hashes, GraphHashData()})
              .first->second);
    return nullptr;
  }

 protected:
  std::map<uint32_t, std::tuple<std::vector<std::vector<int>>,
                                std::vector<uint32_t>, GraphHashData>>
      cache_;
};
