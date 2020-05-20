#include "nauty.hpp"

#include <cmath>

extern "C" {
#include "../third_party/nauty/naugroup.h"
#include "../third_party/nauty/nausparse.h"
#include "../third_party/nauty/nauty.h"
#include "../third_party/nauty/traces.h"
}

Nauty::Nauty(const Graph &G) : G(G) {
  SG_DECL(sg);
  DYNALLSTAT(size_t, sg_v, sg_v_sz);
  DYNALLSTAT(int, sg_d, sg_d_sz);
  DYNALLSTAT(int, sg_e, sg_e_sz);
  DYNALLOC1(size_t, sg_v, sg_v_sz, G.N, "malloc");
  DYNALLOC1(int, sg_d, sg_d_sz, G.N, "malloc");
  DYNALLOC1(int, sg_e, sg_e_sz, 2 * G.M, "malloc");

  sg.nv = G.N;
  sg.nde = 2 * G.M;
  sg.v = sg_v;
  sg.d = sg_d;
  sg.e = sg_e;
  int index_e = 0;
  for (int v = 0; v < G.N; v++) {
    sg.d[v] = G.Adj(v).size();
    sg.v[v] = index_e;
    for (int e : G.Adj(v)) {
      sg.e[index_e++] = e;
    }
  }

  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLOC1(int, lab, lab_sz, G.N, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, G.N, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, G.N, "malloc");

  DEFAULTOPTIONS_SPARSEGRAPH(options);
  statsblk stats;
  // options.getcanon = false;
  options.userautomproc = groupautomproc;
  options.userlevelproc = grouplevelproc;
  sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);

  // Store the group structure.
  group = groupptr(true);

  // DEFAULTOPTIONS_TRACES(options);
  // TracesStats stats;
  // Traces(&sg, lab, ptn, orbits, &options, &stats, NULL);
  num_orbits = stats.numorbits;
  num_generators = stats.numgenerators;
  num_automorphisms = stats.groupsize1 * std::pow(10, stats.groupsize2);

  // Map orbit vertices to orbit indices.
  std::vector<int> orbit_to_index(G.N, G.N);
  int index = 0;
  for (int i = 0; i < G.N; i++) {
    int orbit = orbits[i];
    if (orbit_to_index[orbit] == G.N) orbit_to_index[orbit] = index++;
  }

  // Find a representative with lowest global index for each orbit.
  std::vector<int> orbit_representatives_global(stats.numorbits, full_graph.N);
  orbit_representatives.resize(stats.numorbits, 0);
  this->orbits.resize(G.N, 0);
  for (int i = 0; i < G.N; i++) {
    this->orbits[i] = orbits[i];
    int orbit_index = orbit_to_index[orbits[i]];
    if (G.global[i] < orbit_representatives_global[orbit_index]) {
      orbit_representatives[orbit_index] = i;
      orbit_representatives_global[orbit_index] = G.global[i];
    }
  }
}

Nauty::~Nauty() { freegroup(group); }

std::vector<std::vector<int>> global_automorphisms;
void StoreAutomorhphism(int *p, int n) {
  global_automorphisms.emplace_back();
  global_automorphisms.back().reserve(n);
  for (int i = 0; i < n; i++) global_automorphisms.back().emplace_back(p[i]);
}

const std::vector<std::vector<int>> &Nauty::Automorphisms() {
  if (automorphisms.size()) return automorphisms;
  makecosetreps(group);
  global_automorphisms.reserve(G.N);
  allgroup(group, StoreAutomorhphism);
  automorphisms = std::move(global_automorphisms);
  return automorphisms;
}
