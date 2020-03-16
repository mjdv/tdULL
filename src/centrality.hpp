#pragma once
#include <vector>
#include "graph.hpp"

std::vector<int> DegreeCentrality(const Graph &G);
std::vector<double> BetweennessCentrality(const Graph &G);
std::vector<double> EigenvectorCentrality(const Graph &G, size_t steps = 8);
std::vector<double> PageRankCentrality(const Graph &G, size_t steps = 8, double damping = 0.85);
