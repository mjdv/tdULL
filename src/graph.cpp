#include "graph.hpp"

#include <cassert>
#include <queue>

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

const std::vector<int> &SubGraph::Adj(int v) const {
  assert(v >= 0 && v < vertices.size() && adj.size() == vertices.size());
  return adj[v];
}

std::vector<SubGraph> SubGraph::WithoutVertex(int w) const {
  assert(w >= 0 && w < vertices.size());
  std::vector<SubGraph> cc;
  static std::vector<int> stack;

  // This table will keep the mapping from our indices <-> indices subgraph.
  std::vector<int> sub_vertices;
  sub_vertices.reserve(vertices.size());

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = vertices.size();
  for (int root = 0; root < vertices.size(); ++root) {
    if (vertices_left == 0) break;
    if (root == w) continue;
    if (!vertices[root]->visited) {
      sub_vertices.clear();
      SubGraph component;
      component.vertices.reserve(vertices_left);

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

  static std::vector<int> stack;
  stack.push_back(root);
  vertices[root]->visited = true;
  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();
    result.vertices.emplace_back(vertices[v]);
    for (int nghb : Adj(v))
      if (!vertices[nghb]->visited) {
        stack.push_back(nghb);
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
