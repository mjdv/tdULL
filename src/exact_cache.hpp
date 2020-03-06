#pragma once
#include <array>
#include <map>
#include <vector>

constexpr size_t exactCacheSize = 8;
std::pair<int, int> exactCache(const std::vector<std::vector<int>> &adj);
