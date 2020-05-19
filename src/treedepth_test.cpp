#include "treedepth.hpp"
#include <chrono>
#include <cmath>

#include "centrality.hpp"

int main(int argc, char **argv) {
  // Lower timeout to 5 minutes.
  max_time_treedepth = 5 * 60;
  std::string root = "../input/exact/";
  std::vector<std::pair<std::string, int>> truth_values{
      {"exact_001.gr", 6},  {"exact_003.gr", 11}, {"exact_005.gr", 5},
      {"exact_007.gr", 9},  {"exact_009.gr", 6},  {"exact_011.gr", 5},
      {"exact_013.gr", 7},  {"exact_015.gr", 5},  {"exact_017.gr", 10},
      {"exact_019.gr", 7},  {"exact_021.gr", 5},  {"exact_023.gr", 7},
      {"exact_025.gr", 8},  {"exact_027.gr", 11}, {"exact_029.gr", 12},
      {"exact_031.gr", 11}, {"exact_033.gr", 12}, {"exact_035.gr", 9},
      {"exact_037.gr", 7},  {"exact_039.gr", 7},  {"exact_041.gr", 5},
      {"exact_043.gr", 14}, {"exact_045.gr", 11}, {"exact_047.gr", 25},
      {"exact_049.gr", 10}, {"exact_051.gr", 11}, {"exact_053.gr", 7},
      {"exact_055.gr", 9},  {"exact_057.gr", 13}, {"exact_059.gr", 42},
      {"exact_061.gr", 13}, {"exact_063.gr", 12}, {"exact_065.gr", 11},
      {"exact_067.gr", 9},  {"exact_069.gr", 6},  {"exact_071.gr", 15},
      {"exact_073.gr", 16}, {"exact_077.gr", 12}, {"exact_079.gr", 9},
      {"exact_081.gr", 7},  {"exact_085.gr", 13}, {"exact_087.gr", 13},
      {"exact_089.gr", 15}, {"exact_091.gr", 9},  {"exact_093.gr", 6},
      {"exact_095.gr", 8},  {"exact_097.gr", 7},  {"exact_099.gr", 13},
      {"exact_103.gr", 13}, {"exact_105.gr", 9},  {"exact_107.gr", 11},
      {"exact_109.gr", 9},  {"exact_111.gr", 18}, {"exact_113.gr", 6},
      {"exact_115.gr", 12}, {"exact_125.gr", 10}, {"exact_127.gr", 8},
      {"exact_133.gr", 15}, {"exact_137.gr", 58}, {"exact_141.gr", 13},
      {"exact_145.gr", 7},  {"exact_151.gr", 7},  {"exact_157.gr", 11},
      {"exact_161.gr", 13}, {"exact_165.gr", 9},  {"exact_177.gr", 9},
      {"exact_181.gr", 10}, {"exact_189.gr", 8}};

  auto start_total = std::chrono::steady_clock::now();
  for (auto [fn, true_depth] : truth_values) {
    std::cerr << "Loading example " << fn << "." << std::endl;
    std::ifstream input(root + fn, std::ios::in);
    LoadGraph(input);
    auto start = std::chrono::steady_clock::now();
    int depth = treedepth(full_graph).first;
    if (depth != true_depth) {
      std::cerr << "TEST FAILED!" << std::endl
                << "\t example " << fn << " gives treedepth " << depth
                << " != " << true_depth << std::endl;
      return 1;
    }
    double time_elapsed =
        0.1 * std::round(10 * std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - start)
                                  .count());
    std::cerr << "Example took " << time_elapsed << " seconds." << std::endl
              << std::endl;
    std::cout << fn << "," << depth << "," << time_elapsed << std::endl;
  }
  double time_elapsed =
      0.1 * std::round(10 * std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - start_total)
                                .count());
  std::cerr << "All examples took " << time_elapsed << " seconds." << std::endl;
  std::cout << time_elapsed << std::endl;
}
