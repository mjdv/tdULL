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

  std::istringstream stream_2core(
      "p tdp 11 12 1 2 2 3 3 1 3 4 4 5 5 6 6 7 7 4 6 8 7 9 9 10 10 11");
  LoadGraph(stream_2core);
  auto core = full_graph_as_sub.TwoCore();
  assert(core.vertices.size() == 7);
  assert(core.M == 8);

  std::istringstream stream_3core(
      "p tdp 11 14 1 2 2 3 3 1 3 4 4 5 5 6 6 7 7 4 6 8 7 9 9 10 10 11 4 6 5 7");
  LoadGraph(stream_3core);
  auto cc_core3 = full_graph_as_sub.kCore(3);
  assert(cc_core3.size() == 1);
  auto core3 = cc_core3[0];
  std::cout << core3.vertices.size() << " " << core3.M << std::endl;
  assert(core3.vertices.size() == 4);
  assert(core3.M == 6);

  std::istringstream stream_complement("p tdp 4 4 1 2 2 3 3 4 4 1");
  LoadGraph(stream_complement);
  auto cc_complement = full_graph_as_sub.ComplementComponents();
  std::cout << cc_complement.size() << std::endl;
  for(auto cc : cc_complement)
    std::cout << cc.vertices.size() << std::endl;
  assert(cc_complement.size() == 2);
  assert(cc_complement[0].vertices.size() == 2);
  assert(cc_complement[1].vertices.size() == 2);

  return 0;
}
