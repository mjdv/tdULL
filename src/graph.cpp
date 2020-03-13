#include "graph.hpp"

#include <cassert>

Graph full_graph;
SubGraph full_graph_as_sub;

Graph::Graph(std::istream &stream) {
  std::string str;
  stream >> str;
  assert(str == "p");
  stream >> str;
  assert(str == "tdp");
  stream >> N >> M;

  // Create the vector of vertices.
  vertices.reserve(N);
  for (int v = 0; v < N; v++) vertices.emplace_back(v);
  adj.resize(N);
  for (int e = 0; e < M; e++) {
    int a, b;
    stream >> a >> b;

    // Make them zero indexed.
    a--, b--;

    adj[a].emplace_back(b);
    adj[b].emplace_back(a);
  }
}

SubGraph::SubGraph() : mask(full_graph.N, false) {}

// Create a SubGraph of G with the given (local) vertices
SubGraph::SubGraph(const SubGraph &G, const std::vector<int> &sub_vertices)
    : SubGraph() {
  assert(G.vertices.size() > sub_vertices.size());  // This is silly.
  vertices.reserve(sub_vertices.size());

  // This table will keep the mapping from G indices <-> indices subgraph.
  std::vector<int> new_indices(G.vertices.size(), -1);

  // Add all new vertices to our subgraph.
  for (int v_new = 0; v_new < sub_vertices.size(); ++v_new) {
    int v_old = sub_vertices[v_new];
    new_indices[v_old] = v_new;

    mask[G.vertices[v_old]->n] = true;
    vertices.emplace_back(G.vertices[v_old]);
  }

  // Now find the new adjacency lists.
  for (int v_new = 0; v_new < sub_vertices.size(); ++v_new) {
    int v_old = sub_vertices[v_new];
    std::vector<int> nghbrs;
    nghbrs.reserve(G.Adj(v_old).size());
    for (int nghb_old : G.Adj(v_old))
      if (mask[G.vertices[nghb_old]->n]) {
        assert(new_indices[nghb_old] >= 0 &&
               new_indices[nghb_old] < sub_vertices.size());
        nghbrs.emplace_back(new_indices[nghb_old]);
      }

    // Insert this adjacency list into the subgraph.
    M += nghbrs.size();
    max_degree = std::max(max_degree, nghbrs.size());
    min_degree = std::min(min_degree, nghbrs.size());
    adj.emplace_back(std::move(nghbrs));
  }
  if (min_degree == 0) assert(sub_vertices.size() == 1 && M == 0);
  assert(M % 2 == 0);
  M /= 2;
}

void SubGraph::AssertValidSubGraph() const {
  int N = vertices.size();
  assert(N > 0);
  // First we do some sanity checks on the mask.
  int N_mask = 0;
  for (bool b : mask)
    if (b) N_mask++;
  assert(N_mask == N);
  for (auto vtx : vertices) assert(mask[vtx->n]);

  // Check that no vertex occurs twice.
  assert(std::set<Vertex *>(vertices.begin(), vertices.end()).size() == N);

  // Check that the degrees coincide.
  assert(adj.size() == N);
  size_t max_d = 0, min_d = INT_MAX;
  for (int v = 0; v < N; v++) {
    max_d = std::max(max_d, adj[v].size());
    min_d = std::min(min_d, adj[v].size());
  }
  assert(max_d == max_degree);
  assert(min_d == min_degree);
  assert(max_degree < N);

  // Create a mapping of global to local indices.
  std::map<int, int> glob_2_local;
  for (int v = 0; v < N; ++v) glob_2_local[vertices[v]->n] = v;
  assert(glob_2_local.size() == N);

  // Now check that the adjacency matrix is indeed the `full` adjacency matrix.
  for (int v = 0; v < N; ++v) {
    // First check that there are no doubles in the adj list.
    std::set<int> adj_set(adj[v].begin(), adj[v].end());
    assert(adj_set.size() == adj[v].size());

    // Loop over the global adjacency list.
    int v_glob = vertices[v]->n;
    assert(glob_2_local.count(v_glob));
    for (int nghb_glob : full_graph.adj[v_glob]) {
      assert(glob_2_local.count(nghb_glob) == mask[nghb_glob]);

      // If we contain this nghbour, assert that its also in the local adj
      // list.
      if (mask[nghb_glob]) assert(adj_set.count(glob_2_local.at(nghb_glob)));
    }
  }

  // Now check that this is indeed a connected graph.
  std::vector<int> stack;
  std::vector<bool> visited;
  visited.resize(N, false);
  stack.push_back(0);
  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();
    if (visited[v]) continue;
    visited[v] = true;

    for (int nghb_loc : adj[v])
      if (!visited[nghb_loc]) stack.push_back(nghb_loc);
  }

  for (int v = 0; v < N; ++v) assert(visited[v]);
}

// Gives a vector of all the components of the possibly disconnected subgraph
// of G given by sub_vertices (in local coordinates).
std::vector<SubGraph> SubGraph::ConnectedSubGraphs(
    const std::vector<int> &sub_vertices) const {
  std::vector<bool> in_sub_verts(vertices.size(), false);
  for (int v : sub_vertices) in_sub_verts[v] = true;

  std::vector<SubGraph> cc;

  static std::stack<int> stack;
  static std::vector<int> component;
  component.reserve(sub_vertices.size());
  for (int v : sub_vertices) {
    if (!vertices[v]->visited) {
      component.clear();
      stack.push(v);
      vertices[v]->visited = true;
      while (!stack.empty()) {
        int v = stack.top();
        stack.pop();
        component.push_back(v);
        for (int nghb : Adj(v))
          if (in_sub_verts[nghb] && !vertices[nghb]->visited) {
            stack.push(nghb);
            vertices[nghb]->visited = true;
          }
      }
      cc.emplace_back(*this, component);
    }
  }

  for (int v : sub_vertices) vertices[v]->visited = false;

  return cc;
}

const std::vector<int> &SubGraph::Adj(int v) const {
  assert(v >= 0 && v < vertices.size() && adj.size() == vertices.size());
  return adj[v];
}
int SubGraph::LocalIndex(Vertex *v) const {
  for (int v_local = 0; v_local < vertices.size(); ++v_local) {
    if (vertices[v_local] == v) return v_local;
  }
  assert(false);
}

std::vector<SubGraph> SubGraph::WithoutVertices(
    const std::vector<int> &S) const {
  int N = vertices.size();

  std::vector<bool> in_S(N, false);
  for (auto s : S) in_S[s] = true;

  std::vector<int> remaining;
  remaining.reserve(N - S.size());
  for (int i = 0; i < N; i++) {
    if (!in_S[i]) remaining.push_back(i);
  }

  return ConnectedSubGraphs(remaining);
}

std::vector<SubGraph> SubGraph::WithoutVertex(int w) const {
  assert(w >= 0 && w < vertices.size());
  std::vector<SubGraph> cc;
  static std::vector<int> stack;

  // This table will keep the mapping from our indices <-> indices subgraph.
  static std::vector<int> sub_vertices;
  sub_vertices.reserve(vertices.size());

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = vertices.size();
  for (int root = 0; root < vertices.size(); ++root) {
    if (vertices_left == 0) break;
    if (root == w) continue;
    if (!vertices[root]->visited) {
      sub_vertices.clear();

      // Do a DFS from root, skipping w.
      stack.emplace_back(root);
      vertices[root]->visited = true;
      while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();

        // Add this vertex to the component.
        sub_vertices.push_back(v);

        // Visit all neighbours, skipping w.
        for (int nghb : Adj(v))
          if (!vertices[nghb]->visited && nghb != w) {
            stack.emplace_back(nghb);
            vertices[nghb]->visited = true;
          }
      }

      // Create a SubGraph for this component.
      cc.emplace_back(*this, sub_vertices);
    }
  }
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return cc;
}

std::vector<int> SubGraph::Bfs(int root) const {
  assert(mask[vertices[root]->n]);

  std::vector<int> result;
  result.reserve(vertices.size());

  static std::queue<int> queue;
  queue.push(root);
  vertices[root]->visited = true;
  while (!queue.empty()) {
    int v = queue.front();
    queue.pop();
    result.emplace_back(v);
    for (int nghb : Adj(v))
      if (!vertices[nghb]->visited) {
        queue.push(nghb);
        vertices[nghb]->visited = true;
      }
  }
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return result;
}

SubGraph SubGraph::BfsTree(int root) const {
  assert(mask[vertices[root]->n]);

  SubGraph result;
  result.mask = mask;
  result.vertices.reserve(vertices.size());
  result.adj.resize(vertices.size());

  static std::queue<int> queue;
  queue.push(root);
  vertices[root]->visited = true;
  while (!queue.empty()) {
    int v = queue.front();
    queue.pop();
    result.vertices.emplace_back(vertices[v]);
    for (int nghb : Adj(v))
      if (!vertices[nghb]->visited) {
        queue.push(nghb);
        result.adj[v].push_back(nghb);
        result.adj[nghb].push_back(v);
        vertices[nghb]->visited = true;
      }
  }
  result.M = vertices.size() - 1;
  int total_edges = 0;
  for (int v = 0; v < result.vertices.size(); v++) {
    result.max_degree = std::max(result.max_degree, result.adj[v].size());
    result.min_degree = std::min(result.min_degree, result.adj[v].size());
    total_edges += result.adj[v].size();
  }
  assert(result.min_degree);
  assert(result.M == total_edges / 2);
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return result;
}

SubGraph SubGraph::DfsTree(int root) const {
  assert(mask[vertices[root]->n]);

  SubGraph result;
  result.mask = mask;
  result.vertices.reserve(vertices.size());
  result.adj.resize(vertices.size());

  static std::vector<std::pair<int, int>> stack;
  stack.emplace_back(root, -1);

  while (!stack.empty()) {
    auto [v, prev] = stack.back();
    stack.pop_back();
    if (vertices[v]->visited) continue;
    if (prev != -1) {
      result.adj[v].push_back(prev);
      result.adj[prev].push_back(v);
    }
    vertices[v]->visited = true;
    result.vertices.emplace_back(vertices[v]);
    for (int nghb : Adj(v))
      if (!vertices[nghb]->visited) stack.emplace_back(nghb, v);
  }

  result.M = vertices.size() - 1;
  int total_edges = 0;
  for (int v = 0; v < result.vertices.size(); v++) {
    result.max_degree = std::max(result.max_degree, result.adj[v].size());
    result.min_degree = std::min(result.min_degree, result.adj[v].size());
    total_edges += result.adj[v].size();
  }
  assert(result.min_degree);
  assert(result.M == total_edges / 2);
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return result;
}

std::vector<SubGraph> SubGraph::kCore(int k) const {
  static std::vector<int> stack;
  assert(!IsTreeGraph());
  int N = vertices.size();
  int vertices_left = N;

  // This will keep a list of all the (local) degrees.
  std::vector<int> degrees;
  degrees.reserve(N);
  for (int i = 0; i < N; i++) degrees.push_back(Adj(i).size());

  // Remove all nodes with degree < k.
  for (int v = 0; v < N; v++)
    if (degrees[v] < k && degrees[v] > 0) {
      stack.emplace_back(v);
      degrees[v] = 0;
      while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();
        vertices_left--;

        for (int nghb : Adj(v))
          if (degrees[nghb]) {
            --degrees[nghb];
            if (degrees[nghb] < k && degrees[nghb]) {
              stack.push_back(nghb);
              degrees[nghb] = 0;
            }
          }
      }
    }

  // Nothing was removed, simply return ourself.
  if (vertices_left == N) return {*this};

  // Everything was removed, return the empty graph
  if (vertices_left == 0) return {};

  // Note that the subgraph does not need to be connected, so
  // ceate all subgraphs with vertices that are not removed.
  std::vector<SubGraph> cc;

  static std::vector<int> sub_vertices;
  sub_vertices.reserve(vertices_left);
  for (int v = 0; v < N; v++)
    if (degrees[v] && !vertices[v]->visited) {
      sub_vertices.clear();
      stack.emplace_back(v);
      vertices[v]->visited = true;
      while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();
        sub_vertices.push_back(v);
        for (int nghb : Adj(v))
          if (degrees[nghb] && !vertices[nghb]->visited) {
            stack.push_back(nghb);
            vertices[nghb]->visited = true;
          }
      }
      cc.emplace_back(*this, sub_vertices);
    }

  for (int v = 0; v < N; v++) {
    assert(vertices[v]->visited == (degrees[v] > 0));
    vertices[v]->visited = false;
  }

  return cc;
}

SubGraph SubGraph::TwoCore() const {
  assert(!IsTreeGraph());
  int N = vertices.size();
  int vertices_left = N;

  // This will keep a list of all the (local) degrees.
  std::vector<int> degrees;
  degrees.reserve(N);
  for (int i = 0; i < N; i++) degrees.push_back(Adj(i).size());

  // Remove all leaves, and its adjacent 2 deg nodes.
  for (int v = 0; v < N; v++)
    if (degrees[v] == 1) {
      int cur = v;
      while (degrees[cur] == 1) {
        degrees[cur] = 0;
        vertices_left--;
        for (int nb : Adj(cur)) {
          if (degrees[nb]) {
            degrees[nb]--;
            cur = nb;
            break;
          }
        }
      }
    }

  // Nothing was removed, simply return ourself.
  if (vertices_left == N) return *this;

  // Create subgraph with all vertices that are not removed.
  std::vector<int> sub_vertices;
  sub_vertices.reserve(vertices_left);

  for (int v = 0; v < N; v++)
    if (degrees[v]) sub_vertices.push_back(v);

  return SubGraph(*this, sub_vertices);
}

void LoadGraph(std::istream &stream) {
  full_graph_as_sub = SubGraph();
  full_graph = Graph(stream);

  // Create the subgraph.
  full_graph_as_sub.M = full_graph.M;
  full_graph_as_sub.mask = std::vector<bool>(full_graph.N, true);
  for (int v = 0; v < full_graph.vertices.size(); ++v) {
    Vertex &vertex = full_graph.vertices[v];
    assert(vertex.n == v);
    full_graph_as_sub.vertices.emplace_back(&full_graph.vertices[v]);
    full_graph_as_sub.adj.emplace_back(full_graph.adj[v]);
    full_graph_as_sub.max_degree =
        std::max(full_graph_as_sub.max_degree, full_graph.adj[v].size());
    full_graph_as_sub.min_degree =
        std::min(full_graph_as_sub.min_degree, full_graph.adj[v].size());
  }

  std::cout << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}

SeparatorGenerator::SeparatorGenerator(const SubGraph &G)
    : G(G), N(G.vertices.size()), in_nbh(N, false) {
  // Complete graphs don't have separators. We want this to return a non-empty
  // vector.
  assert(!G.IsCompleteGraph());

  // Datatypes that will be reused.
  std::stack<int> component;
  std::vector<int> separator;
  std::vector<bool> sep_mask(N, false);

  // First we generate all the "seeds": we take the neighborhood of a point
  // (including the point), take all the connected components in the
  // complement, and then take the neighborhoods of those components. Each of
  // those is a minimal separator (in fact, one that separates the original
  // point).
  std::vector<int> neighborhood;
  for (int i = 0; i < N; i++) {
    if (G.Adj(i).size() == N - 1) continue;
    neighborhood.clear();
    neighborhood.push_back(i);
    neighborhood.insert(neighborhood.end(), G.Adj(i).begin(), G.Adj(i).end());

    for (int v : neighborhood) in_nbh[v] = true;

    // Now we start off by enqueueing all neighborhoods of connected components
    // in the complement of this neighborhood.
    for (int j = 0; j < N; j++) {
      if (in_nbh[j] || G.vertices[j]->visited) continue;

      // Reset shared datastructures.
      assert(component.empty());
      separator.clear();

      component.push(j);
      G.vertices[j]->visited = true;

      while (!component.empty()) {
        int cur = component.top();
        component.pop();

        for (int nb : G.Adj(cur)) {
          if (!G.vertices[nb]->visited) {
            if (in_nbh[nb]) {
              separator.push_back(nb);
            } else {
              component.push(nb);
            }
            G.vertices[nb]->visited = true;
          }
        }
      }

      for (auto k : neighborhood) G.vertices[k]->visited = false;

      // std::sort(separator.begin(), separator.end());
      for (int k : separator) sep_mask[k] = true;
      if (done.find(sep_mask) == done.end()) {
        queue.push(separator);
        done.insert(sep_mask);
      }
      for (int k : separator) sep_mask[k] = false;
    }

    for (auto v : G.vertices) v->visited = false;
    for (int v : neighborhood) in_nbh[v] = false;
  }
}

std::vector<Separator> SeparatorGenerator::Next(int k) {
  // Datatypes that will be reused.
  std::stack<int> component;
  std::vector<int> separator;
  std::vector<bool> sep_mask(N, false);

  std::vector<Separator> result;
  while (!queue.empty() && result.size() < k) {
    auto cur_separator = queue.front();
    queue.pop();

    for (int x : cur_separator) {
      for (int j : G.Adj(x)) in_nbh[j] = true;
      for (int j : cur_separator) in_nbh[j] = true;

      for (int j = 0; j < N; j++) {
        if (in_nbh[j] || G.vertices[j]->visited) continue;
        // Reset shared datastructures.
        assert(component.empty());
        separator.clear();

        component.push(j);
        G.vertices[j]->visited = true;

        while (!component.empty()) {
          int cur = component.top();
          component.pop();

          for (int nb : G.Adj(cur)) {
            if (!G.vertices[nb]->visited) {
              if (in_nbh[nb]) {
                separator.push_back(nb);
              } else {
                component.push(nb);
              }
              G.vertices[nb]->visited = true;
            }
          }
        }

        for (auto k : cur_separator) G.vertices[k]->visited = false;
        for (auto k : G.Adj(x)) G.vertices[k]->visited = false;

        // std::sort(separator.begin(), separator.end());
        for (int k : separator) sep_mask[k] = true;
        if (done.find(sep_mask) == done.end()) {
          queue.push(separator);
          done.insert(sep_mask);
        }
        for (int k : separator) sep_mask[k] = false;
      }

      for (int j : cur_separator) in_nbh[j] = false;
      for (int j : G.Adj(x)) in_nbh[j] = false;
      for (int j = 0; j < N; j++) G.vertices[j]->visited = false;
    }

    // Finally, there is the following problem: this algorithm generates all
    // _minimal a,b-separators_, for each pair of vertices a, b in G. These are
    // a strict superset of what we are after: it may be that a set S is a
    // minimal a,b-separator, yet it still has a strict subset which is a
    // separator: it just may be that there is a c,d-separator which is a
    // strict subset.
    //
    // Fortunately, checking whether or not this is the case is relatively
    // fast.

    Separator sep;
    sep.vertices = std::move(cur_separator);
    if (FullyMinimal(sep)) result.push_back(std::move(sep));
  }

  return result;
}

// Checks whether or not the separator of G given by separator is "truly
// minimal": it contains no separator as a strict subset. (Name subject to
// change.)
//
// It also writes some info to the Separator struct: the number of vertices
// and edges in the components that remain when removing this separator from
// the graph.
bool SeparatorGenerator::FullyMinimal(Separator &separator) const {
  // Shared datastructure.
  static std::stack<int> component;

  int N = G.vertices.size();
  std::vector<bool> in_sep(N, false);
  for (int s : separator.vertices) {
    assert(s < N);
    in_sep[s] = true;
  }

  for (int i = 0; i < N; i++) {
    if (!G.vertices[i]->visited && !in_sep[i]) {
      assert(component.empty());

      int comp_N = 1;
      int comp_M = 0;

      component.push(i);
      G.vertices[i]->visited = true;
      while (component.size()) {
        int cur = component.top();
        component.pop();
        for (int nb : G.Adj(cur)) {
          if (!in_sep[nb]) {
            comp_M++;
            if (!G.vertices[nb]->visited) {
              comp_N++;
              component.push(nb);
            }
          }
          G.vertices[nb]->visited = true;
        }
      }
      for (int s : separator.vertices) {
        if (!G.vertices[s]->visited) {
          // This connected component dit not hit this vertex, so we can
          // remove this vertex and have a separator left; so this separator
          // is not fully minimal.
          for (int j = 0; j < N; j++) G.vertices[j]->visited = false;
          return false;
        }
        G.vertices[s]->visited = false;
      }
      separator.comp.push_back({comp_N, comp_M});
    }
  }

  for (int i = 0; i < N; i++) G.vertices[i]->visited = false;

  return true;
}
