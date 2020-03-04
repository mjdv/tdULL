#include "graph.hpp"

#include <cassert>
#include <sstream>

int main() {
  // Load the full graph, this is exact_015.gr.
  std::istringstream stream(
      "p tdp 26 30 18 22 18 13 22 23 22 9 22 2 22 4 22 20 22 24 22 12 10 11 10 "
      "7 10 14 11 1 11 26 11 21 14 26 26 7 2 13 6 17 17 16 17 12 17 19 17 8 17 "
      "15 1 5 1 13 13 21 13 25 13 24 12 3");
  LoadGraph(stream);

  std::cout
      << "Removing one vertex at and check number of connected componenets"
      << std::endl;
  for (int v = 0; v < full_graph_as_sub.vertices.size(); ++v) {
    auto cc = full_graph_as_sub.WithoutVertex(v);

    std::cout << "\tG\\{" << full_graph_as_sub.vertices[v]->n << "} has "
              << cc.size() << " connnected components:" << std::endl;
    for (int c = 0; c < cc.size(); c++)
      std::cout << "\t\tComponent " << c << " has " << cc[c].vertices.size()
                << " vertices and " << cc[c].M << " edges with max_degree "
                << cc[c].max_degree << std::endl;
  }

  return 0;
}
