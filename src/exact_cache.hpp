#pragma once
#include <array>
#include <cassert>
#include <map>
#include <vector>

#define EXACT_CACHE_SIZE 6
constexpr size_t exactCacheSize = EXACT_CACHE_SIZE;
std::pair<int, int> exactCache(const std::vector<std::vector<int>> &adj);
