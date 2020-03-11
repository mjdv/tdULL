#include "treedepth.hpp"

#include "centrality.hpp"

int main(int argc, char **argv) {
  std::string root = "../input/exact/";
  std::vector<std::pair<std::string, int>> truth_values{
      {"exact_001.gr", 6},  {"exact_003.gr", 11}, {"exact_005.gr", 5},
      {"exact_007.gr", 9},  {"exact_009.gr", 6},  {"exact_011.gr", 5},
      {"exact_013.gr", 7},  {"exact_015.gr", 5},  {"exact_017.gr", 10},
      {"exact_019.gr", 7},  {"exact_021.gr", 5},  {"exact_023.gr", 7},
      {"exact_025.gr", 8},  {"exact_027.gr", 11}, {"exact_029.gr", 12},
      {"exact_031.gr", 11}, {"exact_033.gr", 12}, {"exact_035.gr", 9},
      {"exact_037.gr", 7},  {"exact_039.gr", 7},  {"exact_041.gr", 5},
      {"exact_045.gr", 11}, {"exact_049.gr", 10}, {"exact_053.gr", 7},
      {"exact_055.gr", 9},  {"exact_067.gr", 9},  {"exact_079.gr", 9},
      {"exact_081.gr", 7},  {"exact_091.gr", 9},  {"exact_093.gr", 6},
      {"exact_095.gr", 8},  {"exact_097.gr", 7},  {"exact_105.gr", 9},
      {"exact_109.gr", 9},  {"exact_113.gr", 6},  {"exact_125.gr", 10},
      {"exact_127.gr", 8},  {"exact_145.gr", 7},  {"exact_151.gr", 7},
      {"exact_165.gr", 9},  {"exact_177.gr", 9},  {"exact_181.gr", 10}};

  time_t start_total, end_total;
  time(&start_total);
  for (auto [fn, true_depth] : truth_values) {
    time_t start, end;
    std::cout << "Loading example " << fn << "." << std::endl;
    time(&start);
    std::ifstream input(root + fn, std::ios::in);
    LoadGraph(input);
    auto [depth, tree] = treedepth(full_graph_as_sub);
    if (depth != true_depth) {
      std::cout << "TEST FAILED!" << std::endl
                << "\t example " << fn << " gives treedepth " << depth
                << " != " << true_depth << std::endl;
      return 1;
    }
    time(&end);
    std::cout << "Example took " << difftime(end, start) << " seconds."
              << std::endl
              << std::endl;
  }
  time(&end_total);
  std::cout << "All examples took " << difftime(end_total, start_total)
            << " seconds." << std::endl;
}
