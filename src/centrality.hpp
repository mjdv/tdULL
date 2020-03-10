#pragma once
#include <vector>
#include "graph.hpp"

std::vector<int> DegreeCentrality(const SubGraph &G);
std::vector<double> BetweennessCentrality(const SubGraph &G);
