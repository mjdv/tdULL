#include "treedepth.hpp"

int main(int argc, char **argv) {
  std::string root = "../input/exact/";
  std::vector<std::pair<std::string, int>> truth_values{
      {"exact_001.gr", 6}, {"exact_003.gr", 11}, {"exact_005.gr", 5},
      {"exact_007.gr", 9}, {"exact_009.gr", 6},  {"exact_011.gr", 5},
      {"exact_013.gr", 7}, {"exact_015.gr", 5},  {"exact_017.gr", 10},
      {"exact_019.gr", 7}, {"exact_021.gr", 5},  {"exact_023.gr", 7},
      {"exact_025.gr", 8}, {"exact_027.gr", 11}, {"exact_029.gr", 12},
      {"exact_031.gr", 11}};

  time_t start_total, end_total;
  time(&start_total);
  for (auto [fn, true_depth] : truth_values) {
    time_t start, end;
    time(&start);
    std::ifstream input(root + fn, std::ios::in);
    LoadGraph(input);
    int depth = treedepth(full_graph_as_sub);
    if (depth != true_depth) {
      std::cout << "TEST FAILED!" << std::endl
                << "\t example " << fn << " gives treedepth " << depth
                << " != " << true_depth << std::endl;
      return 1;
    }
    time(&end);
    std::cout << "Example " << fn << " took " << difftime(end, start)
              << " seconds." << std::endl;
  }
  time(&end_total);
  std::cout << "All examples took " << difftime(end_total, start_total)
            << " seconds." << std::endl;
}
