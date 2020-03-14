#pragma once
#include <vector>
#include "graph.hpp"

std::vector<int> DegreeCentrality(const SubGraph &G);
std::vector<double> BetweennessCentrality(const SubGraph &G);
std::vector<double> EigenvectorCentrality(const SubGraph &G, size_t steps = 8);
std::vector<double> PageRankCentrality(const SubGraph &G, size_t steps = 8, double damping = 0.85);
