#pragma once
#include <cstdint>
#include <vector>

// Calculates the graph hash for graph G given by its adj lists.
// Evaluates a maximum of N steps of hash updating.
std::pair<uint32_t, std::vector<uint32_t>> GraphHash(
    const std::vector<std::vector<int>> &G);

// Heuristically calculates a isomoprhism between two graphs on basis of hashes.
std::pair<bool, std::vector<int>> GraphHashIsomphismMapping(
    const std::vector<std::vector<int>> &G1,
    const std::vector<std::vector<int>> &G2);

// Calculates whether two graphs are an isomoprhism based on the hashes.
bool GraphIsomorphism(const std::vector<std::vector<int>> &G1,
                      const std::vector<std::vector<int>> &G2,
                      const std::vector<int> &mapping);

// Uses a graph hash as heuristic for checking if two graphs are isomporhm.
bool GraphHashIsomorphism(const std::vector<std::vector<int>> &G1,
                          const std::vector<std::vector<int>> &G2);
