#include "nauty.hpp"

#include <sstream>

int main() {
  // Load a complete graph.
  std::istringstream complete_6(
      "p tdp 6 15 1 2 1 3 1 4 1 5 1 6 2 3 2 4 2 5 2 6 3 4 3 5 3 6 4 5 4 6 5 6");
  LoadGraph(complete_6);

  Nauty nauty_kn_6(full_graph);
  nauty_kn_6.Automorphisms();
  assert(nauty_kn_6.num_orbits == 1);
  assert(nauty_kn_6.num_automorphisms == 720);
  assert(nauty_kn_6.orbit_representatives.size() == 1);
  assert(nauty_kn_6.automorphisms.size() == nauty_kn_6.num_automorphisms);

  // Load a line graph
  std::istringstream line_6("p tdp 6 5 1 2 2 3 3 4 4 5 5 6");
  LoadGraph(line_6);
  Nauty nauty_line_6(full_graph);
  nauty_line_6.Automorphisms();
  assert(nauty_line_6.num_orbits == 3);
  assert(nauty_line_6.num_automorphisms == 2);
  assert(nauty_line_6.orbit_representatives.size() == nauty_line_6.num_orbits);
  assert(nauty_line_6.automorphisms.size() == nauty_line_6.num_automorphisms);

  return 0;
}
