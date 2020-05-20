#include "separator.hpp"

#include <cassert>

// Initializes a separator of G. This checks whether or not the separator
// of G given by the vertices is "truly minimal": it contains no separator
// as a strict subset. (Name subject to change.)
//
// It also writes some info to the Separator struct: the number of vertices
// and edges in the components that remain when removing this separator from
// the graph.
Separator::Separator(const Graph &G, const std::vector<int> &vertices)
    : vertices(vertices), fully_minimal(true) {
  // Shared datastructure.
  static std::stack<int> component;

  std::vector<bool> visited(G.N, false);
  std::vector<bool> in_sep(G.N, false);
  for (int s : vertices) {
    assert(s < G.N);
    in_sep[s] = true;
  }

  for (int i = 0; i < G.N; i++) {
    if (!visited[i] && !in_sep[i]) {
      assert(component.empty());

      int comp_N = 1;
      int comp_M = 0;

      component.push(i);
      visited[i] = true;
      while (component.size()) {
        int cur = component.top();
        component.pop();
        for (int nb : G.Adj(cur)) {
          if (!in_sep[nb]) {
            comp_M++;
            if (!visited[nb]) {
              comp_N++;
              component.push(nb);
            }
          }
          visited[nb] = true;
        }
      }
      for (int s : vertices) {
        if (!visited[s]) {
          // This connected component dit not hit this vertex, so we can
          // remove this vertex and have a separator left; so this separator
          // is not fully minimal.
          fully_minimal = false;
        }
        visited[s] = false;
      }
      // components.emplace_back(comp_N, comp_M);
      largest_component = std::max(largest_component, {comp_N, comp_M / 2});
    }
  }
}

SeparatorGenerator::SeparatorGenerator(const Graph &G)
    : G(G), in_nbh(G.N, false), sep_mask(G.N, false) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;
  static std::vector<int> neighborhood;

  // Complete graphs don't have separators. We want this to return a
  // non-empty vector.
  assert(!G.IsCompleteGraph());
  if (G.N > 30) {
    nauty = std::unique_ptr<Nauty>(new Nauty(G));
    if (nauty->num_automorphisms > 1 && nauty->num_automorphisms < 1000000) {
      const auto &automorphisms = nauty->Automorphisms();
      in_automorphisms.resize(G.N);
      for (int i = 0; i < automorphisms.size(); i++) {
        const auto &f = automorphisms[i];
        for (int v = 0; v < f.size(); v++)
          if (v != f[v]) in_automorphisms[v].emplace_back(i);
      }
    } else
      nauty.reset();
  }

  // First we generate all the "seeds": we take the neighborhood of a
  // point (including the point), take all the connected components in
  // the complement, and then take the neighborhoods of those
  // components. Each of those is a minimal separator (in fact, one that
  // separates the original point).
  std::vector<bool> visited(G.N, false);
  for (int i = 0; i < G.N; i++) {
    if (G.Adj(i).size() == G.N - 1) continue;
    neighborhood.clear();
    neighborhood.push_back(i);
    neighborhood.insert(neighborhood.end(), G.Adj(i).begin(), G.Adj(i).end());

    for (int v : neighborhood) in_nbh[v] = true;

    // Now we start off by enqueueing all neighborhoods of connected
    // components in the complement of this neighborhood.
    for (int j = 0; j < G.N; j++) {
      if (in_nbh[j] || visited[j]) continue;

      // Reset shared datastructures.
      assert(component.empty());
      separator.clear();

      component.push(j);
      visited[j] = true;

      while (!component.empty()) {
        int cur = component.top();
        component.pop();

        for (int nb : G.Adj(cur)) {
          if (!visited[nb]) {
            if (in_nbh[nb]) {
              separator.push_back(nb);
            } else {
              component.push(nb);
            }
            visited[nb] = true;
          }
        }
      }

      for (auto k : neighborhood) visited[k] = false;

      // std::sort(separator.begin(), separator.end());
      for (int k : separator) sep_mask[k] = true;
      if (done.find(sep_mask) == done.end()) {
        queue.push(separator);
        done.insert(sep_mask);

        Separator sep(G, separator);
        if (sep.fully_minimal) buffer.emplace_back(std::move(sep));
      }
      for (int k : separator) sep_mask[k] = false;
    }

    for (int j = 0; j < G.N; j++) visited[j] = false;
    for (int v : neighborhood) in_nbh[v] = false;
  }
}

std::vector<Separator> SeparatorGenerator::Next(int k) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;

  std::vector<bool> visited(G.N, false);
  std::vector<bool> sep_mask_automorphism;
  std::unordered_set<int> try_automorphisms;
  if (nauty) sep_mask_automorphism = std::vector<bool>(G.N, false);

  while (!queue.empty() && buffer.size() < k) {
    auto cur_separator = queue.front();
    queue.pop();

    for (int x : cur_separator) {
      for (int j : G.Adj(x)) in_nbh[j] = true;
      for (int j : cur_separator) in_nbh[j] = true;

      for (int j = 0; j < G.N; j++) {
        if (in_nbh[j] || visited[j]) continue;
        // Reset shared datastructures.
        assert(component.empty());
        separator.clear();

        component.push(j);
        visited[j] = true;

        while (!component.empty()) {
          int cur = component.top();
          component.pop();

          for (int nb : G.Adj(cur)) {
            if (!visited[nb]) {
              if (in_nbh[nb]) {
                separator.push_back(nb);
              } else {
                component.push(nb);
              }
              visited[nb] = true;
            }
          }
        }

        for (auto k : cur_separator) visited[k] = false;
        for (auto k : G.Adj(x)) visited[k] = false;

        // Check whether we already have the separator.
        bool sep_found = false;
        for (int k : separator) sep_mask[k] = true;
        sep_found = done.find(sep_mask) != done.end();

        if (nauty && !sep_found) {
          // If we have automorphisms, check these as well.
          try_automorphisms.clear();
          for (int k : separator)
            try_automorphisms.insert(in_automorphisms[k].begin(),
                                     in_automorphisms[k].end());

          for (int i : try_automorphisms) {
            const auto &fn = nauty->automorphisms[i];
            for (int k : separator) sep_mask_automorphism[fn[k]] = true;
            sep_found = done.find(sep_mask_automorphism) != done.end();
            for (int k : separator) sep_mask_automorphism[fn[k]] = false;
            if (sep_found) break;
          }
        }

        // If we have not found the separator, add it.
        if (!sep_found) {
          queue.push(separator);
          done.insert(sep_mask);

          Separator sep(G, separator);
          if (sep.fully_minimal) buffer.emplace_back(std::move(sep));
        }
        for (int k : separator) sep_mask[k] = false;
      }

      for (int j : cur_separator) in_nbh[j] = false;
      for (int j : G.Adj(x)) in_nbh[j] = false;
      for (int j = 0; j < G.N; j++) visited[j] = false;
    }
  }

  return std::move(buffer);
}
