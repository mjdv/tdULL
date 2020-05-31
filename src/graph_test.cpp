#include "graph.hpp"

#include <cassert>
#include <sstream>

#include "separator.hpp"

// Helper function to verify that the given seperators are indeed fully minimal,
// and that the sizes of the components are correct :-).
void TestSeparators(const Graph &G, const std::vector<Separator> &seps) {
  for (const auto &sep : seps) {
    // First check that the given graph is indeed fully minimal.
    Graph H = G;
    for (int i = 1; i < sep.vertices.size(); i++) {
      int v = sep.vertices[i - 1];
      // Get the subgraph after removing vertex v.
      auto cc = H.WithoutVertex(H.LocalIndex(G.global[v]));
      assert(cc.size() == 1);
      H = cc[0];
    }

    // Remove all vertices and assert that this does decompose the graph.
    auto cc = G.WithoutVertices(sep.vertices);
    assert(cc.size() > 1);

    // Put the sizes of these components into a vector.
    std::vector<std::pair<int, int>> cc_sizes;
    for (auto H : cc) cc_sizes.emplace_back(H.N, H.M);

    // Find the largest component
    auto cc_largest = *std::max_element(cc_sizes.begin(), cc_sizes.end());

    // compare to the separator sizes
    assert(cc_largest == sep.largest_component);
  }
}

// Helper function to verify that the given articulationpoints function is
// correct.
void TestArticulationPoints(const Graph &G) {
  std::vector<int> aps_real;
  // Create a poors mans version of articulat point checker.
  for (int v = 0; v < G.N; v++)
    if (G.WithoutVertex(v).size() > 1) aps_real.push_back(v);

  std::sort(aps_real.begin(), aps_real.end());

  std::vector<int> aps = G.ArticulationPoints();
  std::sort(aps.begin(), aps.end());

  assert(aps == aps_real);
}

void TestSymmetricNeighboorhoods(const Graph &G_original) {
  auto [G_new, vertices_original] = G_original.WithoutSymmetricNeighboorhoods();
  G_new.AssertValidGraph();
  assert(vertices_original.size() == G_new.N);
  // Nothing happened..
  if (vertices_original.size() == G_original.N)
    for (int v = 0; v < G_original.N; v++) {
      assert(vertices_original[v].size() == 1);
      assert(vertices_original[v][0] == v);
    }

  // Complete graphs are an exception.
  if (G_new.IsCompleteGraph()) return;

  // Check that there are no double neighboorhoods in G_new.
  for (int v = 0; v < G_new.N; v++) {
    for (int w = v + 1; w < G_new.N; w++) {
      std::set<int> adj_w(G_new.Adj(w).begin(), G_new.Adj(w).end());
      std::set<int> adj_v(G_new.Adj(v).begin(), G_new.Adj(v).end());
      adj_w.erase(v);
      adj_v.erase(w);
      assert(adj_v != adj_w);
    }
  }
}

int main() {
  // Load the full graph, this is exact_015.gr.
  std::istringstream stream_015(
      "p tdp 26 30 18 22 18 13 22 23 22 9 22 2 22 4 22 20 22 24 22 12 10 11 "
      "10 "
      "7 10 14 11 1 11 26 11 21 14 26 26 7 2 13 6 17 17 16 17 12 17 19 17 8 "
      "17 "
      "15 1 5 1 13 13 21 13 25 13 24 12 3");
  LoadGraph(stream_015);

  //  Test Graph.WithoutVertex.
  std::cout
      << "Removing one vertex at and check number of connected componenets"
      << std::endl;
  for (int v = 0; v < full_graph.N; ++v) {
    auto cc = full_graph.WithoutVertex(v);

    std::cout << "\tG\\{" << full_graph.global[v] << "} has " << cc.size()
              << " connnected components:" << std::endl;
    for (int c = 0; c < cc.size(); c++)
      std::cout << "\t\tComponent " << c << " has " << cc[c].N
                << " vertices and " << cc[c].M << " edges with max_degree "
                << cc[c].max_degree << std::endl;
  }
  auto gen15 = SeparatorGenerator(full_graph);
  auto v_ams15 = gen15.Next(1'000'000);
  assert(v_ams15.size() == 8);
  TestSeparators(full_graph, v_ams15);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  std::istringstream stream_2core(
      "p tdp 11 12 1 2 2 3 3 1 3 4 4 5 5 6 6 7 7 4 6 8 7 9 9 10 10 11");
  LoadGraph(stream_2core);
  auto core = full_graph.TwoCore();
  assert(core.N == 7);
  assert(core.M == 8);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  std::istringstream stream_3core(
      "p tdp 11 14 1 2 2 3 3 1 3 4 4 5 5 6 6 7 7 4 6 8 7 9 9 10 10 11 4 6 5 "
      "7");
  LoadGraph(stream_3core);
  auto cc_core3 = full_graph.kCore(3);
  assert(cc_core3.size() == 1);
  auto core3 = cc_core3[0];
  std::cout << core3.N << " " << core3.M << std::endl;
  assert(core3.N == 4);
  assert(core3.M == 6);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  std::istringstream stream_allminsep("p tdp 6 5 1 2 1 3 1 4 1 5 1 6");
  LoadGraph(stream_allminsep);
  auto gen = SeparatorGenerator(full_graph);
  auto v_ams = gen.Next(1'000'000);
  TestSeparators(full_graph, v_ams);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  assert(v_ams.size() == 1);
  assert(v_ams[0].vertices.size() == 1);

  std::istringstream stream_allminsep2("p tdp 6 6 1 2 2 3 3 4 4 5 5 6 6 1");
  LoadGraph(stream_allminsep2);
  auto gen2 = SeparatorGenerator(full_graph);
  auto v_ams2 = gen2.Next(1'000'000);
  assert(v_ams2.size() == 9);
  TestSeparators(full_graph, v_ams2);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  std::cout << "The minimal separators of the 6-cycle are:" << std::endl;
  for (auto v : v_ams2) {
    assert(v.vertices.size() == 2);
    for (int x : v.vertices) std::cout << full_graph.global[x] << " ";
    std::cout << std::endl;
  }

  std::istringstream stream_ams3("p tdp 6 6 1 2 2 3 3 4 4 1 3 5 4 6");
  LoadGraph(stream_ams3);
  auto gen3 = SeparatorGenerator(full_graph);
  auto v_ams3 = gen3.Next(1'000'000);
  TestSeparators(full_graph, v_ams3);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);

  std::cout << "The minimal separators of the 4-cycle with two extra leaves "
               "attached to adjacent nodes are:"
            << std::endl;
  for (auto v : v_ams3) {
    for (int x : v.vertices) std::cout << full_graph.global[x] << " ";
    std::cout << std::endl;
  }

  std::cout << "Contracting one edge of the four-cycle...\n";
  Graph my_graph = full_graph.Contract({0, 1});
  std::cout << "This gives graph:\n";
  for (int v = 0; v < my_graph.N; v++) {
    std::cout << my_graph.global[v] << "(" << v << "): ";
    for (auto nb : my_graph.Adj(v)) std::cout << my_graph.global[nb] << " ";
    std::cout << std::endl;
  }

  std::cout << "Contracting one more edge ...\n";
  my_graph = my_graph.Contract({0, 1});
  std::cout << "This gives graph:\n";
  for (int v = 0; v < my_graph.N; v++) {
    std::cout << my_graph.global[v] << "(" << v << "): ";
    for (auto nb : my_graph.Adj(v)) std::cout << my_graph.global[nb] << " ";
    std::cout << std::endl;
  }

  // Graph exact_043.gr.
  std::istringstream stream_043(
      "p tdp 40 129 9 26 9 21 9 14 9 19 9 8 9 25 9 29 9 11 9 10 9 7 26 11 26 "
      "14 "
      "26 10 26 8 26 27 24 30 24 14 24 38 30 18 30 20 30 5 30 1 30 2 30 37 "
      "30 "
      "23 "
      "30 28 30 3 14 17 14 15 14 25 14 23 14 27 14 13 17 13 17 16 5 31 5 2 5 "
      "4 "
      "5 "
      "1 5 3 5 13 5 29 5 22 31 1 31 37 31 18 31 39 31 2 7 12 7 25 7 19 7 27 "
      "12 "
      "34 12 1 12 2 12 13 12 33 6 28 6 1 6 2 6 13 6 3 28 13 33 36 33 39 33 "
      "38 "
      "33 "
      "37 33 13 33 35 36 29 25 23 25 22 25 8 25 10 25 11 25 29 25 19 25 21 "
      "25 "
      "40 "
      "13 32 13 16 13 29 13 8 13 4 13 2 13 15 13 3 13 1 2 22 2 32 2 29 2 3 2 "
      "4 "
      "22 3 22 18 22 40 22 38 22 1 22 21 34 35 34 37 34 38 34 39 35 32 1 3 1 "
      "32 "
      "1 4 1 29 18 20 18 29 32 38 32 16 32 39 32 37 19 29 19 20 29 20 29 3 "
      "21 "
      "38 "
      "20 38 20 40 38 23 38 27 11 27 11 8 27 8 27 10 23 40 8 10");
  LoadGraph(stream_043);
  std::cout << "Computing all minimal separators for exact_043.gr...\n";
  auto gen43 = SeparatorGenerator(full_graph);
  auto v_ams43 = gen43.Next(1'000'000);
  assert(v_ams43.size() == 664);
  TestSeparators(full_graph, v_ams43);
  TestArticulationPoints(full_graph);
  TestSymmetricNeighboorhoods(full_graph);
  std::cout << "There are " << v_ams43.size() << " of them." << std::endl;

  return 0;
}
