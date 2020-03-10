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

  std::istringstream stream_allminsep(
      "p tdp 6 5 1 2 1 3 1 4 1 5 1 6");
  LoadGraph(stream_allminsep);
  auto v_ams = full_graph_as_sub.AllMinimalSeparators();
  assert(v_ams.size() == 1);
  assert(v_ams[0].size() == 1);
  assert(full_graph_as_sub.vertices[v_ams[0][0]]->n == 0);

  std::istringstream stream_allminsep2(
      "p tdp 6 6 1 2 2 3 3 4 4 5 5 6 6 1");
  LoadGraph(stream_allminsep2);
  auto v_ams2 = full_graph_as_sub.AllMinimalSeparators();
  assert(v_ams2.size() == 9);

  std::cout << "The minimal separators of the 6-cycle are:" << std::endl;
  for(auto v : v_ams2) {
    assert(v.size() == 2);
    for(int x : v)
      std::cout << full_graph_as_sub.vertices[x]->n << " ";
    std::cout << std::endl;
  }

  // Graph exact_043.gr.
  std::istringstream stream_043(
      "p tdp 40 129 9 26 9 21 9 14 9 19 9 8 9 25 9 29 9 11 9 10 9 7 26 11 26 14 "
      "26 10 26 8 26 27 24 30 24 14 24 38 30 18 30 20 30 5 30 1 30 2 30 37 30 23 "
      "30 28 30 3 14 17 14 15 14 25 14 23 14 27 14 13 17 13 17 16 5 31 5 2 5 4 5 "
      "1 5 3 5 13 5 29 5 22 31 1 31 37 31 18 31 39 31 2 7 12 7 25 7 19 7 27 12 "
      "34 12 1 12 2 12 13 12 33 6 28 6 1 6 2 6 13 6 3 28 13 33 36 33 39 33 38 33 "
      "37 33 13 33 35 36 29 25 23 25 22 25 8 25 10 25 11 25 29 25 19 25 21 25 40 "
      "13 32 13 16 13 29 13 8 13 4 13 2 13 15 13 3 13 1 2 22 2 32 2 29 2 3 2 4 "
      "22 3 22 18 22 40 22 38 22 1 22 21 34 35 34 37 34 38 34 39 35 32 1 3 1 32 "
      "1 4 1 29 18 20 18 29 32 38 32 16 32 39 32 37 19 29 19 20 29 20 29 3 21 38 "
      "20 38 20 40 38 23 38 27 11 27 11 8 27 8 27 10 23 40 8 10");
  LoadGraph(stream_043);
  std::cout << "Computing all minimal separators for exact_043.gr...\n";
  auto v_ams43 = full_graph_as_sub.AllMinimalSeparators();
  std::cout << "There are " << v_ams43.size() << " of them." << std::endl;
  return 0;
}
