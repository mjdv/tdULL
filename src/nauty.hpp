#include "graph.hpp"

struct group_struct;
struct Nauty {
  Nauty(const Graph &G);
  ~Nauty();
  const std::vector<std::vector<int>> &Automorphisms();

  // Stats.
  size_t num_automorphisms = 0;
  size_t num_generators = 0;
  size_t num_orbits = 0;

  // Calculated values.
  std::vector<int> orbit_representatives;
  std::vector<std::vector<int>> automorphisms;

 protected:
  const Graph &G;
  group_struct *group = nullptr;
};
