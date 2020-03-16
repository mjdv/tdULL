#include "graph.hpp"
#include "centrality.hpp"

#include <cassert>
#include <sstream>

int main() {

  
  // Create graph K_{2, 3}.
  std::istringstream bipartite_stream(
      "p tdp 5 6\n"
      "1 4\n"
      "1 5\n"
      "2 4\n"
      "2 5\n"
      "3 4\n"
      "3 5\n");
  LoadGraph(bipartite_stream);

  std::cout
      << "Degree centrality of K_{3, 2}:"
      << std::endl;
  auto degree = DegreeCentrality(full_graph);
  for (int v = 0; v < full_graph.N; ++v) {
      std::cout << degree[v] << " ";
  }
  std::cout << std::endl;

  std::cout
      << "Betweenness centrality of K_{3, 2}:"
      << std::endl;
  auto betweenness = BetweennessCentrality(full_graph);
  for (int v = 0; v < full_graph.N; ++v) {
      std::cout << betweenness[v] << " ";
  }
  std::cout << std::endl;

  std::cout << std::endl;

  // Create graph L_7
  std::istringstream line_stream(
      "p tdp 7 6\n"
      "1 2\n"
      "2 3\n"
      "3 4\n"
      "4 5\n"
      "5 6\n"
      "6 7\n");
  LoadGraph(line_stream);

  std::cout
      << "Degree centrality of L_7:"
      << std::endl;
  degree = DegreeCentrality(full_graph);
  for (int v = 0; v < full_graph.N; ++v) {
      std::cout << degree[v] << " ";
  }
  std::cout << std::endl;

  std::cout
      << "Betweenness centrality of L_7:"
      << std::endl;
  betweenness = BetweennessCentrality(full_graph);
  for (int v = 0; v < full_graph.N; ++v) {
      std::cout << betweenness[v] << " ";
  }
  std::cout << std::endl;

  
  return 0;
}
