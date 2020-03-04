#include "treedepth.hpp"

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Expecting 2 arguments." << std::endl;
    return 1;
  }

  std::ifstream input;
  input.open(argv[1], std::ios::in);
  LoadGraph(input);
  input.close();

  time_t start, end;
  time(&start);
  std::cout << "treedepth gives result:" << std::endl;
  std::cout << treedepth(full_graph_as_sub) << std::endl;
  time(&end);
  std::cout << "Elapsed time is " << difftime(end, start) << " seconds.\n";
}
