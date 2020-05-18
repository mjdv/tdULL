#include <chrono>
#include <cmath>
#include <exception>

#include "treedepth.hpp"

std::string ExtractFileName(const std::string& str) {
  return str.substr(str.find_last_of("/") + 1);
}

std::vector<int> NonDominatedOrbits(const Graph& G, const Nauty& nauty) {
  std::vector<int> result;
  result.reserve(nauty.num_orbits);
  std::vector<bool> in_nbh(G.N, false);
  std::vector<bool> dominated(G.N, false);
  for (int v_prime = 0; v_prime < G.N; v_prime++) {
    bool contract = false;
    for (int nb : G.adj[v_prime]) in_nbh[nb] = true;
    for (int v : nauty.orbit_representatives) {
      if (G.adj[v].size() > G.adj[v_prime].size()) continue;
      // We want v' < v.
      if (G.adj[v].size() == G.adj[v_prime].size() && v_prime >= v) continue;
      if (nauty.orbits[v] == nauty.orbits[v_prime]) continue;
      if (dominated[v]) continue;
      bool subset = true;
      for (int nb : G.adj[v])
        if (!in_nbh[nb] && nb != v_prime) {
          subset = false;
          break;
        }
      if (subset) {
        // N(v) \ v' \subset N(v') \ v
        dominated[v] = true;
      }
    }
    for (int nb : G.adj[v_prime]) in_nbh[nb] = false;
  }
  for (int v : nauty.orbit_representatives)
    if (!dominated[v]) result.emplace_back(v);
  return result;
}

int main(int argc, char** argv) {
  LoadGraph(std::cin);
  Nauty nauty_full(full_graph);
  Nauty nauty_full_contract(full_graph.WithoutSymmetricNeighboorhoods().first);
  auto non_dominated_orbits = NonDominatedOrbits(full_graph, nauty_full);
  size_t leaves = 0;
  for (int v : non_dominated_orbits)
    if (full_graph.Adj(v).size() == 1) leaves++;
  std::cerr << nauty_full.num_orbits << "," << non_dominated_orbits.size()
            << "," << non_dominated_orbits.size() - leaves << std::endl;
  //            << nauty_full.num_orbits << ","
  //            << nauty_full_contract.num_automorphisms << std::endl;
  return 0;

  auto start = std::chrono::steady_clock::now();
  try {
    //    auto seperator_gen = SeparatorGenerator(full_graph_as_sub);
    //    size_t total_count = 0;
    //    while (seperator_gen.HasNext()) {
    //      total_count += seperator_gen.Next(100000).size();
    //      time(&end);
    //      std::cerr << "Total number of separators is " << total_count
    //                << ". Speed is " << double(total_count) / difftime(end,
    //                start)
    //                << " seps / s.\n";
    //    }

    auto [td, tree] = treedepth(full_graph);
    double time_elapsed =
        0.1 * std::round(10 * std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - start)
                                  .count());
    std::cerr << "Treedepth is: " << td << std::endl;
    std::cerr << "Elapsed time is " << time_elapsed << " seconds.\n";

    std::cout << td << std::endl;
    for (int parent : tree) std::cout << parent << std::endl;

    std::cerr << td << "," << time_elapsed << ", " << std::endl;
    return 0;
  } catch (std::exception& e) {
    double time_elapsed =
        0.1 * std::round(10 * std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - start)
                                  .count());
    std::cerr << "Failed! Encountered exception:\"" << e.what() << "\"."
              << std::endl;
    std::cerr << -1 << "," << time_elapsed << "," << e.what() << std::endl;
    return 1;
  }
}
