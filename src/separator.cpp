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
    : G(G), in_nbh(G.N, false), sep_mask(G.N, false) {}

SeparatorGeneratorDirected::SeparatorGeneratorDirected(const Graph &G,
                                                       const int source)
    : SeparatorGenerator(G), source(source) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;
  static std::vector<int> neighborhood;

  // If the chosen source is adjacent to everything, we necessarily want it as
  // root, as there won't be separators without it.
  if (G.Adj(source).size() == G.N - 1) return;

  // First we generate all the "seeds": we take the neighborhood of the
  // source, take all the connected components in the complement, and then
  // take the neighborhoods of those components. Each of those is a minimal
  // separator (in fact, one that separates the source).
  std::vector<bool> visited(G.N, false);

  neighborhood.clear();
  neighborhood.push_back(source);
  neighborhood.insert(neighborhood.end(), G.Adj(source).begin(),
                      G.Adj(source).end());

  for (int v : neighborhood) in_nbh[v] = true;

  // Now we start off by enqueueing all neighborhoods of connected
  // components in the complement of this neighborhood.
  for (int j = 0; j < G.N; j++) {
    // First we find the component of in_nbh that we want.
    if (in_nbh[j] || visited[j]) continue;

    // Reset shared datastructures.
    assert(component.empty());
    separator.clear();

    component.push(j);
    visited[j] = true;

    std::vector<bool> in_separator(G.N, false);

    while (!component.empty()) {
      int cur = component.top();
      component.pop();

      for (int nb : G.Adj(cur)) {
        if (!visited[nb]) {
          if (in_nbh[nb]) {
            in_separator[nb] = true;
            separator.push_back(nb);
          } else {
            component.push(nb);
          }
          visited[nb] = true;
        }
      }
    }

    // Then we check what the component of the source in the complement is.
    std::vector<bool> source_comp(G.N, false);
    source_comp[source] = true;
    component.push(source);

    while (!component.empty()) {
      int cur = component.top();
      component.pop();

      for (int nb : G.Adj(cur)) {
        if (!source_comp[nb] && !in_separator[nb]) {
          source_comp[nb] = true;
          component.push(nb);
        }
      }
    }

    for (auto k : neighborhood) visited[k] = false;

    // std::sort(separator.begin(), separator.end());
    for (int k : separator) sep_mask[k] = true;
    if (done.find(sep_mask) == done.end()) {
      queue.push(make_pair(separator, source_comp));
      done.insert(sep_mask);

      Separator sep(G, separator);
      if (sep.fully_minimal) {
        buffer.emplace_back(std::move(sep));
      }
    }
    for (int k : separator) sep_mask[k] = false;
  }
}

std::vector<Separator> SeparatorGeneratorDirected::Next(int k) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;

  std::vector<bool> visited(G.N, false);

  while (!queue.empty() && buffer.size() < k) {
    auto [cur_separator, cur_source_comp] = queue.front();
    queue.pop();

    for (int x : cur_separator) {
      for (int j : cur_separator) in_nbh[j] = true;

      // If we can only reach source_comp vertices from here, skip.
      bool any_change = false;
      for (int j : G.Adj(x)) {
        if (!in_nbh[j] && !cur_source_comp[j]) {
          in_nbh[j] = true;
          any_change = true;
        }
      }
      if (!any_change) {
        for (int j : cur_separator) in_nbh[j] = false;
        continue;
      }

      // We are only going to start a dfs from a node if it is not in
      // cur_source_comp. (Might be able to push this just a tiny bit faster
      // by restricting j to neighbors of the nodes added to nbh.)
      for (int j = 0; j < G.N; j++) {
        if (in_nbh[j] || visited[j] || cur_source_comp[j]) continue;
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

        for (int k : separator) sep_mask[k] = true;
        if (done.find(sep_mask) == done.end()) {
          std::vector<bool> new_source_comp(cur_source_comp);
          for (int removed : cur_separator) {
            if (sep_mask[removed]) continue;
            for (int nb : G.Adj(removed)) {
              if (cur_source_comp[nb]) {
                component.push(removed);
                new_source_comp[removed] = true;
                break;
              }
            }
          }
          while (!component.empty()) {
            int cur = component.top();
            component.pop();
            for (int nb : G.Adj(cur)) {
              if (!new_source_comp[nb] && !sep_mask[nb]) {
                new_source_comp[nb] = true;
                component.push(nb);
              }
            }
          }

          queue.push(make_pair(separator, new_source_comp));
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

SeparatorGeneratorUndirected::SeparatorGeneratorUndirected(const Graph &G)
    : SeparatorGenerator(G) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;
  static std::vector<int> neighborhood;

  // Complete graphs don't have separators. We want this to return a
  // non-empty vector.


  // Temporarily removing this assert for 2Core testing.
  // assert(!G.IsCompleteGraph());
  if(G.IsCompleteGraph())
    return;

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

std::vector<Separator> SeparatorGeneratorUndirected::Next(int k) {
  // Datatypes that will be reused.
  static std::stack<int> component;
  static std::vector<int> separator;

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

SeparatorGenerator2Core::SeparatorGenerator2Core(const Graph &G) :
  SeparatorGenerator(G), two_core(G.kCore(2)[0]), two_core_gen(two_core) {

  std::map<int, int> G_to_local;
  std::map<int, int> two_core_to_local;
  G_to_two_core = std::vector<int>(G.N, -1);
  two_core_to_G = std::vector<int>(two_core.N, -1);

  for(int i = 0; i < G.N; i++) {
    G_to_local[G.global[i]] = i;
  }
  for(int i = 0; i < two_core.N; i++) {
    two_core_to_local[two_core.global[i]] = i;
  }

  for(int i = 0; i < two_core.N; i++) {
    two_core_to_G[i] = G_to_local[two_core.global[i]];
    G_to_two_core[G_to_local[two_core.global[i]]] = i;
  }

  std::vector<bool> in_two_core(G.N, false);
  two_core_adj = std::vector<bool>(G.N, false);
  for(int i = 0; i < G.N; i++) {
    if(G_to_two_core[i] != -1) {
      in_two_core[i] = true;
      two_core_adj[i] = true;
      for(int nb : G.Adj(i))
        two_core_adj[nb] = true;
    }
  }

  for(int i = 0; i < G.N; i++) {
    if(two_core_adj[i]) {
      buffer.push_back(Separator(G, {i}));
    }
  }
}

bool SeparatorGenerator2Core::HasNext() {
  if(!buffer.empty())
    return true;
  if(!two_core_gen.HasNext())
    return false;
  auto separators = Next(1);
  if(separators.empty())
    return false;
  buffer.insert(buffer.end(), separators.begin(), separators.end());
  return true;
}

std::vector<Separator> SeparatorGenerator2Core::Next(int k) {
  while(buffer.size() < k && two_core_gen.HasNext()) {
    auto two_core_separators = two_core_gen.Next(2*k);
    for(auto sep : two_core_separators) {
      bool contains_two_core_adj = false;
      for(auto v : sep.vertices) {
        if(two_core_adj[two_core_to_G[v]]) {
          contains_two_core_adj = true;
          break;
        }
      }
      if(contains_two_core_adj)
        continue;

      std::vector<int> sep_in_G_coords(sep.vertices.size());
      for(int i = 0; i < sep.vertices.size(); i++)
        sep_in_G_coords[i] = two_core_to_G[sep.vertices[i]];
      buffer.emplace_back(G, sep_in_G_coords);
    }
  }
  return std::move(buffer);
}
