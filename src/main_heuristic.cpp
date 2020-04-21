#include <exception>
#include <chrono>
#include <cmath>

#include "treedepth.hpp"

std::string ExtractFileName(const std::string& str) {
  return str.substr(str.find_last_of("/") + 1);
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Expecting 1 argument." << std::endl;
    return 1;
  }

  std::ifstream input;
  input.open(argv[1], std::ios::in);
  LoadGraph(input);
  input.close();

  treedepth_heuristic(full_graph);

  return 0;
}
