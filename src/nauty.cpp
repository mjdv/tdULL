#include "nauty.hpp"

extern "C" {
#include "../third_party/nauty/naugroup.h"
#include "../third_party/nauty/nausparse.h"
#include "../third_party/nauty/nauty.h"
#include "../third_party/nauty/traces.h"
}

// void automproc(int count, int *perm, int *orbits, int numorbits, int
// stabvertex,
//               int n) {
//  std::cout << count << " " << n << std::endl;
//}

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

  orbit_representatives.resize(stats.numorbits, 0);
  int index = 0;
  for (int i = 0; i < G.N; i++) {
    if (orbits[i] == i) orbit_representatives[index++] = i;
  }
}

Nauty::~Nauty() { freegroup(group); }

std::vector<std::vector<int>> global_automorphisms;
void StoreAutomorhphism(int *p, int n) {
  global_automorphisms.emplace_back();
  global_automorphisms.back().assign(n, *p);
}

const std::vector<std::vector<int>> &Nauty::Automorphisms() {
  if (automorphisms.size()) return automorphisms;
  makecosetreps(group);
  global_automorphisms.reserve(G.N);
  allgroup(group, StoreAutomorhphism);
  automorphisms = std::move(global_automorphisms);
  return automorphisms;
}
