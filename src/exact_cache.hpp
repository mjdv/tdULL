#pragma once
#include <array>
#include <cassert>
#include <map>
#include <vector>

constexpr size_t exactCacheSize = 0;
std::pair<int, int> exactCache(const std::vector<std::vector<int>> &adj);
