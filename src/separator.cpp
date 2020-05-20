#include "separator.hpp"

#include <cassert>

// Initializes a separator of G. This checks whether or not the separator
// of G given by the vertices is "truly minimal": it contains no separator
// as a strict subset. (Name subject to change.)
//
// It also writes some info to the Separator struct: the number of vertices
// and edges in the components that remain when removing this separator from
// the graph.
Separator SeparatorGenerator::InitSeparator(const std::vector<int> &vertices) {
  // Shared datastructures.
  static std::stack<int> component;

  Separator result;
  std::vector<bool> visited(G.N, false);
  std::vector<bool> in_sep(G.N, false);

  // We must expand the vertices.
  result.vertices.reserve(vertices.size());
  for (int s : vertices) {
    assert(s < G.N);
    in_sep[s] = true;
    result.vertices.insert(result.vertices.end(), vertices_original[s].begin(),
                           vertices_original[s].end());
  }

  assert(vertices_original.size() == G.N);
  for (int i = 0; i < G.N; i++) {
    if (!visited[i] && !in_sep[i]) {
      assert(component.empty());

      int comp_N = vertices_original[i].size();
      int comp_M = 0;

      component.push(i);
      visited[i] = true;
      while (component.size()) {
        int cur = component.top();
        component.pop();
        for (int nb : G.Adj(cur)) {
          if (!in_sep[nb]) {
            comp_M +=
                vertices_original[cur].size() * vertices_original[nb].size();
            if (!visited[nb]) {
              comp_N += vertices_original[nb].size();
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
          result.fully_minimal = false;
        }
        visited[s] = false;
      }
      // components.emplace_back(comp_N, comp_M);
      result.largest_component =
          std::max(result.largest_component, {comp_N, comp_M / 2});
    }
  }
  return result;
}

SeparatorGenerator::SeparatorGenerator(const Graph &G_original)
    : G_original(G_original),
      in_nbh(G_original.N, false),
      sep_mask(G_original.N, false) {
  // Complete graphs don't have separators. We want this to return a
  // non-empty vector.
  assert(!G_original.IsCompleteGraph());

  // First we contract the graph, the old indices are stored inside the
  // vertices_original member variable.
  std::tie(G, vertices_original) = G.WithoutSymmetricNeighboorhoods();
  if (G_original.N == full_graph.N)
    std::cerr << "full_graph: contracted graph has " << G.N << " / "
              << G_original.N << " vertices." << std::endl;

  // Minor edge case.
  if (G.IsCompleteGraph()) {
    G = G_original;
    vertices_original.resize(G.N);
    for (int v = 0; v < G.N; v++) vertices_original[v] = {v};
  }

  // Next we generate all the "seeds": we take the neighborhood of a
  // point (including the point), take all the connected components in
  // the complement, and then take the neighborhoods of those
  // components. Each of those is a minimal separator (in fact, one that
  // separates the original point).
  {
    // Datatypes that will be reused.
    static std::stack<int> component;
    static std::vector<int> separator;
    static std::vector<int> neighborhood;
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

          auto sep = InitSeparator(separator);
          if (sep.fully_minimal) buffer.emplace_back(std::move(sep));
        }
        for (int k : separator) sep_mask[k] = false;
      }

      for (int j = 0; j < G.N; j++) visited[j] = false;
      for (int v : neighborhood) in_nbh[v] = false;
    }
  }
}

std::vector<Separator> SeparatorGenerator::Next(int k) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;
  static std::vector<int> vertices_expanded;

  std::vector<bool> visited(G.N, false);

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

        // std::sort(separator.begin(), separator.end());
        for (int k : separator) sep_mask[k] = true;
        if (done.find(sep_mask) == done.end()) {
          queue.push(separator);
          done.insert(sep_mask);

          auto sep = InitSeparator(separator);
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
