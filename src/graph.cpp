#include "graph.hpp"

#include <cassert>

Graph full_graph;
std::vector<bool> full_graph_mask;
std::vector<std::vector<int>> global_to_vertices;
std::map<std::vector<int>, int> vertices_to_global;

int GetIndex(std::vector<int> &vertices);

Graph::Graph(std::istream &stream) {
  std::string str;
  stream >> str;
  assert(str == "p");
  stream >> str;
  assert(str == "tdp");
  stream >> N >> M;

  // Create the vector of vertices.
  global.reserve(N);
  for (int v = 0; v < N; v++) global.push_back(v);
  adj.resize(N);
  for (int e = 0; e < M; e++) {
    int a, b;
    stream >> a >> b;

    // Make them zero indexed.
    a--, b--;

    adj[a].emplace_back(b);
    adj[b].emplace_back(a);
  }

  for (int v = 0; v < N; v++) {
    min_degree = std::min(min_degree, adj[v].size());
    max_degree = std::max(max_degree, adj[v].size());
  }
  full_graph_mask.resize(N, false);
}

Graph::Graph() {}

// Create a Graph of G with the given (local) vertices
Graph::Graph(const Graph &G, std::vector<int> &sub_vertices) : Graph() {
  assert(G.N > sub_vertices.size());  // This is silly.
  N = sub_vertices.size();
  global.reserve(sub_vertices.size());
  adj.reserve(sub_vertices.size());

  if (G.N > full_graph_mask.size()) full_graph_mask.resize(G.N, false);
  for (int v : sub_vertices) full_graph_mask[G.global[v]] = true;

  // This table will keep the mapping from G indices <-> indices subgraph.
  std::vector<int> new_indices(G.N, -1);

  // Add all new vertices to our subgraph, in order.
  for (int v_new = 0; v_new < sub_vertices.size(); ++v_new) {
    int v_old = sub_vertices[v_new];
    new_indices[v_old] = v_new;

    global.emplace_back(G.global[v_old]);
  }

  // Now find the new adjacency lists.
  for (int v_new = 0; v_new < sub_vertices.size(); ++v_new) {
    int v_old = sub_vertices[v_new];
    std::vector<int> nghbrs;
    nghbrs.reserve(G.Adj(v_old).size());
    for (int nghb_old : G.Adj(v_old))
      if (full_graph_mask[G.global[nghb_old]]) {
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

  for (int v : sub_vertices) full_graph_mask[G.global[v]] = false;
}

void Graph::AssertValidGraph() const {
  assert(N > 0);

  // Check that no vertex occurs twice.
  assert(std::set<int>(global.begin(), global.end()).size() == N);

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
  for (int v = 0; v < N; ++v) glob_2_local[global[v]] = v;
  assert(glob_2_local.size() == N);

  // Now check that the adjacency matrix is indeed the `full` adjacency matrix.
  for (int v = 0; v < N; ++v) {
    // First check that there are no doubles in the adj list.
    std::set<int> adj_set(adj[v].begin(), adj[v].end());
    assert(adj_set.size() == adj[v].size());

    // Loop over the global adjacency list.
    int v_glob = global[v];
    assert(glob_2_local.count(v_glob));
    /*for (int nghb_glob : full_graph.adj[v_glob]) {
      assert(glob_2_local.count(nghb_glob) == mask[nghb_glob]);

      // If we contain this nghbour, assert that its also in the local adj
      // list.
      if (mask[nghb_glob]) assert(adj_set.count(glob_2_local.at(nghb_glob)));
    }*/
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
std::vector<Graph> Graph::ConnectedGraphs(
    const std::vector<int> &sub_vertices) const {
  std::vector<bool> in_sub_verts(N, false);
  std::vector<bool> visited(N, false);
  for (int v : sub_vertices) in_sub_verts[v] = true;

  std::vector<Graph> cc;

  static std::stack<int> stack;
  static std::vector<int> component;
  component.reserve(sub_vertices.size());
  for (int v : sub_vertices) {
    if (!visited[v]) {
      component.clear();
      stack.push(v);
      visited[v] = true;
      while (!stack.empty()) {
        int v = stack.top();
        stack.pop();
        component.push_back(v);
        for (int nghb : Adj(v))
          if (in_sub_verts[nghb] && !visited[nghb]) {
            stack.push(nghb);
            visited[nghb] = true;
          }
      }
      cc.emplace_back(*this, component);
    }
  }

  return cc;
}

const std::vector<int> &Graph::Adj(int v) const {
  assert(v >= 0 && v < N && adj.size() == N);
  return adj[v];
}

int Graph::LocalIndex(int global_index) const {
  for (int v_local = 0; v_local < N; ++v_local) {
    if (global[v_local] == global_index) return v_local;
  }
  assert(false);
}

bool Graph::ConnectedSubset(const std::vector<int> vertices) const {
  static std::stack<int> s;
  assert(s.empty());
  assert(vertices.size());

  std::vector<bool> in_vertices(N, false);
  for (int v : vertices) in_vertices[v] = true;

  std::vector<bool> visited(N, false);
  s.push(vertices[0]);
  visited[vertices[0]] = true;
  while (s.size()) {
    int cur = s.top();
    s.pop();
    for (int nb : Adj(cur)) {
      if (!visited[nb] && in_vertices[nb]) {
        visited[nb] = true;
        s.push(nb);
      }
    }
  }
  for (int v : vertices) {
    if (!visited[v]) return false;
  }
  return true;
}

Graph Graph::Contract(const std::vector<int> &contractors) const {
  assert(ConnectedSubset(contractors));

  std::vector<bool> in_contractors(N, false);
  for (auto v : contractors) in_contractors[v] = true;

  Graph contracted;
  contracted.N = N - ((int)contractors.size()) + 1;
  contracted.global.reserve(contracted.N);

  // Mappings between contracted and old indices.
  std::vector<int> local_to_contracted_local(N, -1);
  std::vector<int> contracted_local_to_local(contracted.N, -1);

  for (int i = 0; i < N; i++) {
    if (!in_contractors[i]) {
      local_to_contracted_local[i] = contracted.global.size();
      contracted_local_to_local[contracted.global.size()] = i;

      contracted.global.push_back(global[i]);
    }
  }

  std::vector<int> original_contractors;
  for (int contractor : contractors)
    for (int original_contractor : global_to_vertices[global[contractor]])
      original_contractors.push_back(original_contractor);

  contracted.global.push_back(GetIndex(original_contractors));
  contracted.adj.resize(contracted.N);

  // Adjacency list for non-contracted vertices.
  for (int i = 0; i < contracted.N - 1; i++) {
    bool connected_to_contractors = false;
    contracted.adj[i].reserve(Adj(contracted_local_to_local[i]).size());
    for (int nb : Adj(contracted_local_to_local[i])) {
      if (in_contractors[nb]) {
        if (!connected_to_contractors) {
          connected_to_contractors = true;
          contracted.adj[i].push_back(contracted.N - 1);
          contracted.adj[contracted.N - 1].push_back(i);
          contracted.M += 2;
        }
      } else {
        contracted.adj[i].push_back(local_to_contracted_local[nb]);
        contracted.M++;
      }
    }
  }

  assert(contracted.M % 2 == 0);
  contracted.M /= 2;

  for (int i = 0; i < contracted.N; i++) {
    contracted.min_degree =
        std::min(contracted.min_degree, contracted.adj[i].size());
    contracted.max_degree =
        std::max(contracted.max_degree, contracted.adj[i].size());
  }

  return contracted;
}

std::vector<Graph> Graph::WithoutVertices(const std::vector<int> &S) const {
  std::vector<bool> in_S(N, false);
  for (auto s : S) in_S[s] = true;

  std::vector<int> remaining;
  remaining.reserve(N - S.size());
  for (int i = 0; i < N; i++) {
    if (!in_S[i]) remaining.push_back(i);
  }

  return ConnectedGraphs(remaining);
}

std::vector<Graph> Graph::WithoutVertex(int w) const {
  assert(w >= 0 && w < N);
  std::vector<Graph> cc;
  static std::vector<int> stack;

  // This table will keep the mapping from our indices <-> indices subgraph.
  static std::vector<int> sub_vertices;
  sub_vertices.reserve(N);
  std::vector<bool> visited(N, false);

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = N;
  for (int root = 0; root < N; ++root) {
    if (vertices_left == 0) break;
    if (root == w) continue;
    if (!visited[root]) {
      sub_vertices.clear();

      // Do a DFS from root, skipping w.
      stack.emplace_back(root);
      visited[root] = true;
      while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();

        // Add this vertex to the component.
        sub_vertices.push_back(v);

        // Visit all neighbours, skipping w.
        for (int nghb : Adj(v))
          if (!visited[nghb] && nghb != w) {
            stack.emplace_back(nghb);
            visited[nghb] = true;
          }
      }

      // Create a Graph for this component.
      cc.emplace_back(*this, sub_vertices);
    }
  }
  // Reset the visited field.
  return cc;
}

std::pair<Graph, std::vector<std::vector<int>>>
Graph::WithoutSymmetricNeighboorhoods() const {
  // Find all pairs v, w with N(v) = N(w) or N(v)\{w} = N(w)\{v}.
  std::vector<int> vertices_contract;
  std::vector<std::vector<int>> vertices_original;
  std::vector<bool> contracted(N, false);
  std::vector<bool> in_nbh(N, false);

  vertices_contract.reserve(N);
  vertices_original.reserve(N);
  for (int v = 0; v < N; v++) {
    if (contracted[v]) continue;
    vertices_contract.emplace_back(v);
    vertices_original.emplace_back(std::vector<int>{v});
    for (int nb : adj[v]) in_nbh[nb] = true;
    for (int w = v + 1; w < N; w++) {
      if (adj[v].size() != adj[w].size() || contracted[w]) continue;

      // Check if the neighbours of w coincide with that of v.
      bool contract = true;
      for (int nb : adj[w])
        if (!in_nbh[nb] && nb != v) {
          contract = false;
          break;
        }
      if (contract) {
        vertices_original.back().emplace_back(w);
        contracted[w] = true;
      }
    }
    for (int nb : adj[v]) in_nbh[nb] = false;
  }

  if (vertices_contract.size() < N)
    return {Graph(*this, vertices_contract), std::move(vertices_original)};
  else
    return {*this, std::move(vertices_original)};
}

std::vector<int> Graph::NonDominatedVertices(
    const std::vector<int> &vertices) const {
  std::vector<int> result;
  result.reserve(vertices.size());
  std::vector<bool> in_nbh(N, false);
  std::vector<bool> dominated(N, false);
  for (int v_prime : vertices) {
    bool contract = false;
    for (int nb : adj[v_prime]) in_nbh[nb] = true;
    for (int v : vertices) {
      if (adj[v].size() > adj[v_prime].size()) continue;
      // We want v' < v.
      if (adj[v].size() == adj[v_prime].size() && v_prime >= v) continue;
      if (dominated[v]) continue;
      bool subset = true;
      for (int nb : adj[v])
        if (!in_nbh[nb] && nb != v_prime) {
          subset = false;
          break;
        }
      if (subset) {
        // N(v) \ v' \subset N(v') \ v
        dominated[v] = true;
      }
    }
    for (int nb : adj[v_prime]) in_nbh[nb] = false;
  }
  for (int v : vertices)
    if (!dominated[v]) result.emplace_back(v);
  return result;
}

std::vector<int> Graph::Bfs(int root) const {
  // assert(mask[vertices[root]->n]);
  std::vector<bool> visited(N, false);

  std::vector<int> result;
  result.reserve(N);

  static std::queue<int> queue;
  queue.push(root);
  visited[root] = true;
  while (!queue.empty()) {
    int v = queue.front();
    queue.pop();
    result.emplace_back(v);
    for (int nghb : Adj(v))
      if (!visited[nghb]) {
        queue.push(nghb);
        visited[nghb] = true;
      }
  }
  return result;
}

Graph Graph::BfsTree(int root) const {
  // assert(mask[vertices[root]->n]);

  std::vector<bool> visited(N, false);

  Graph result;
  result.N = N;
  // result.mask = mask;
  result.global.reserve(N);
  result.adj.resize(N);

  static std::queue<int> queue;
  queue.push(root);
  visited[root] = true;
  while (!queue.empty()) {
    int v = queue.front();
    queue.pop();
    result.global.emplace_back(global[v]);
    for (int nghb : Adj(v))
      if (!visited[nghb]) {
        queue.push(nghb);
        result.adj[v].push_back(nghb);
        result.adj[nghb].push_back(v);
        visited[nghb] = true;
      }
  }
  result.M = N - 1;
  int total_edges = 0;
  for (int v = 0; v < result.N; v++) {
    result.max_degree = std::max(result.max_degree, result.adj[v].size());
    result.min_degree = std::min(result.min_degree, result.adj[v].size());
    total_edges += result.adj[v].size();
  }
  assert(result.min_degree);
  assert(result.M == total_edges / 2);
  return result;
}

Graph Graph::DfsTree(int root) const {
  // assert(mask[vertices[root]->n]);
  std::vector<bool> visited(N, false);

  Graph result;
  result.N = N;
  result.global.reserve(N);
  result.adj.resize(N);

  static std::vector<std::pair<int, int>> stack;
  stack.emplace_back(root, -1);

  while (!stack.empty()) {
    auto [v, prev] = stack.back();
    stack.pop_back();
    if (visited[v]) continue;
    if (prev != -1) {
      result.adj[v].push_back(prev);
      result.adj[prev].push_back(v);
    }
    visited[v] = true;
    result.global.emplace_back(global[v]);
    for (int nghb : Adj(v))
      if (!visited[nghb]) stack.emplace_back(nghb, v);
  }

  result.M = N - 1;
  int total_edges = 0;
  for (int v = 0; v < result.N; v++) {
    result.max_degree = std::max(result.max_degree, result.adj[v].size());
    result.min_degree = std::min(result.min_degree, result.adj[v].size());
    total_edges += result.adj[v].size();
  }
  assert(result.min_degree);
  assert(result.M == total_edges / 2);
  return result;
}

void ArticulationPointsHelper(const Graph &G, int u, int &d,
                              std::vector<bool> &visited,
                              std::vector<int> &depth, std::vector<int> &low,
                              std::vector<int> &parent,
                              std::vector<int> &result) {
  int children = 0;
  visited[u] = true;

  depth[u] = d;
  low[u] = d;
  d++;

  bool is_ap = false;

  for (int v : G.Adj(u)) {
    if (!visited[v]) {
      children++;
      parent[v] = u;
      ArticulationPointsHelper(G, v, d, visited, depth, low, parent, result);
      low[u] = std::min(low[u], low[v]);
      if (low[v] >= depth[u]) is_ap = true;
    } else if (v != parent[u]) {
      low[u] = std::min(low[u], depth[v]);
    }
  }
  if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && is_ap))
    result.push_back(u);
}

std::vector<int> Graph::ArticulationPoints() const {
  std::vector<int> parent(N, -1);
  std::vector<int> depth(N, INT_MAX);
  std::vector<int> low(N, INT_MAX);
  std::vector<bool> visited(N, false);
  int d = 0;

  std::vector<int> result;
  for (int u = 0; u < N; u++) {
    if (!visited[u])
      ArticulationPointsHelper(*this, u, d, visited, depth, low, parent,
                               result);
  }

  return result;
}

// void writeautom(int *p, int n) {
//  for (int i = 0; i < n; i++) {
//    fprintf(stderr, " %2d", p[i]);
//  }
//  fprintf(stderr, "\n");
//}

// std::vector<int> Graph::nauty_call() const {
//  // std::cerr << "start nauty call" << std::endl;
//  sparsegraph g = this->sparsegraph_nauty();
//
//  // std::cerr << "done with sparsegraph" << std::endl;
//
//  int *lab = (int *)calloc(this->N, sizeof(int));
//  int *ptn = (int *)calloc(this->N, sizeof(int));
//  int *orbits = (int *)calloc(this->N, sizeof(int));
//
//  // grouprec *group;
//
//  DEFAULTOPTIONS_SPARSEGRAPH(options);
//  options.writeautoms = false;
//  options.writemarkers = false;
//  options.outfile = stderr;
//  options.userautomproc = groupautomproc;
//  options.userlevelproc = grouplevelproc;
//
//  statsblk stats;
//  sparsenauty(&g, lab, ptn, orbits, &options, &stats, NULL);
//
//  /*
//  std::cerr << "done with nauty" << std::endl;
//  std::cerr << "\nResults" << std::endl;
//  std::cerr << "groupsize : ";
//  writegroupsize(stderr, stats.grpsize1, stats.grpsize2);
//  std::cerr << std::endl;
//
//  std::cerr << "nr of orbits : " << stats.numorbits << std::endl;
//  std::cerr << "Nr of generators : " << stats.numgenerators << std::endl;
//
//  std::cerr << "ORBITS" << std::endl;
//  for(int i=0; i<this->N; i++){
//    int cur_orbit = *(orbits + i);
//    std::cerr << cur_orbit << std::endl;
//  }
//  std::cerr << "END OF ORBITS" << std::endl;
//  std::cerr << "All group elements" << std::endl;
//  */
//
//  std::vector<int> orbit_representatives = std::vector<int>(stats.numorbits,
//  0); int index = 0; for (int i = 0; i < this->N; i++) {
//    if (orbits[i] == i) orbit_representatives[index++] = i;
//  }
//
//  // group = groupptr(false);
//  // makecosetreps(group);
//  // allgroup(group, writeautom);
//
//  // std::cerr << "end of nauty_call" << std::endl;
//
//  return orbit_representatives;
//}

std::vector<Graph> Graph::kCore(int k) const {
  static std::vector<int> stack;
  assert(!IsTreeGraph());
  int vertices_left = N;

  // This will keep a list of all the (local) degrees.
  static std::vector<int> degrees;
  if (degrees.size() < N) degrees.resize(N);
  for (int i = 0; i < N; i++) degrees[i] = Adj(i).size();

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
  std::vector<int> remaining;
  remaining.reserve(vertices_left);
  for (int v = 0; v < N; v++)
    if (degrees[v]) remaining.emplace_back(v);

  return ConnectedGraphs(remaining);
}

Graph Graph::TwoCore() const {
  assert(!IsTreeGraph());
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

  return Graph(*this, sub_vertices);
}

void LoadGraph(std::istream &stream) {
  full_graph = Graph(stream);

  global_to_vertices = std::vector<std::vector<int>>();
  vertices_to_global = std::map<std::vector<int>, int>();

  for (int i = 0; i < full_graph.N; i++) {
    std::vector<int> i_vec = {i};
    global_to_vertices.push_back(i_vec);
    vertices_to_global[i_vec] = i;
  }

  std::cerr << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}

int GetIndex(std::vector<int> &vertices) {
  sort(vertices.begin(), vertices.end());
  auto it = vertices_to_global.find(vertices);
  if (it == vertices_to_global.end()) {
    int new_index = vertices_to_global.size();
    vertices_to_global[vertices] = new_index;
    global_to_vertices.push_back(vertices);
    full_graph_mask.push_back(false);
    return new_index;
  }
  return it->second;
}
