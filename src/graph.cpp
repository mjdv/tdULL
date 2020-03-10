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

const std::vector<int> &SubGraph::Adj(int v) const {
  assert(v >= 0 && v < vertices.size() && adj.size() == vertices.size());
  return adj[v];
}

std::vector<SubGraph> SubGraph::WithoutVertex(int w) const {
  assert(w >= 0 && w < vertices.size());
  std::vector<SubGraph> cc;
  static std::vector<int> stack;

  // This table will keep the mapping from our indices <-> indices subgraph.
  std::vector<int> new_indices(vertices.size(), -1);
  std::vector<int> old_indices(vertices.size(), -1);

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = vertices.size();
  for (int root = 0; root < vertices.size(); ++root) {
    if (vertices_left == 0) break;
    if (root == w) continue;
    if (!vertices[root]->visited) {
      SubGraph component;
      component.vertices.reserve(vertices_left);

      // Do a DFS from root, skipping w.
      stack.emplace_back(root);
      vertices[root]->visited = true;
      while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();

        // Create index mapping between current graph and subgraph.
        new_indices[v] = component.vertices.size();
        old_indices[component.vertices.size()] = v;

        // Insert this vertex into the component.
        component.vertices.emplace_back(vertices[v]);
        component.mask[vertices[v]->n] = true;
        vertices_left--;

        // Visit all neighbours, skipping w.
        for (int nghb : Adj(v))
          if (!vertices[nghb]->visited && nghb != w) {
            stack.emplace_back(nghb);
            vertices[nghb]->visited = true;
          }
      }

      // Find the (trimmed) adjacency list for each vertex.
      component.adj.reserve(component.vertices.size());
      for (int v_new = 0; v_new < component.vertices.size(); ++v_new) {
        int v_old = old_indices[v_new];
        assert(v_old >= 0 && v_old < vertices.size());

        // Find the trimmed adjacency list for this vertex.
        std::vector<int> nghbrs;
        nghbrs.reserve(Adj(v_old).size());
        for (int nghb : Adj(v_old))
          if (nghb != w) {
            assert(new_indices[nghb] >= 0 &&
                   new_indices[nghb] < component.vertices.size());
            nghbrs.emplace_back(new_indices[nghb]);
          }

        // Insert this adjacency list into the component.
        component.M += nghbrs.size();
        component.max_degree = std::max(component.max_degree, nghbrs.size());
        component.adj.emplace_back(std::move(nghbrs));
      }

      assert(component.M % 2 == 0);
      component.M /= 2;
      cc.emplace_back(std::move(component));
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
    total_edges += result.adj[v].size();
  }
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
    total_edges += result.adj[v].size();
  }
  assert(result.M == total_edges / 2);
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return result;
}

SubGraph SubGraph::TwoCore() const {
  assert(!IsTreeGraph());
  int N = vertices.size();

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
        for (int nb : Adj(cur)) {
          if (degrees[nb]) {
            degrees[nb]--;
            cur = nb;
            break;
          }
        }
      }
    }

  // Create subgraph with all vertices that are not removed.
  SubGraph H;
  H.vertices.reserve(N);

  // This table will keep the mapping from our indices <-> indices subgraph.
  std::vector<int> new_indices(vertices.size(), -1);
  std::vector<int> old_indices(vertices.size(), -1);

  for (int v = 0; v < N; v++) {
    assert(degrees[v] >= 0 && degrees[v] != 1);
    if (degrees[v]) {
      // Create index mapping between current graph and subgraph.
      new_indices[v] = H.vertices.size();
      old_indices[H.vertices.size()] = v;

      // Add the vertex to the subgraph.
      H.vertices.emplace_back(vertices[v]);
      H.mask[vertices[v]->n] = true;
    }
  }
  assert(H.vertices.size());

  // Now find the trimmed adjacency lists.
  for (int v_old = 0; v_old < N; v_old++) {
    if (degrees[v_old]) {
      std::vector<int> nghbrs;
      nghbrs.reserve(Adj(v_old).size());
      for (int nghb : Adj(v_old))
        if (degrees[nghb]) {
          assert(new_indices[nghb] >= 0 &&
                 new_indices[nghb] < H.vertices.size());
          nghbrs.emplace_back(new_indices[nghb]);
        }

      // Insert this adjacency list into the subgraph.
      H.M += nghbrs.size();
      H.max_degree = std::max(H.max_degree, nghbrs.size());
      H.adj.emplace_back(std::move(nghbrs));
    }
  }
  assert(H.vertices.size() < vertices.size());
  assert(H.M % 2 == 0);
  H.M /= 2;
  return H;
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
  }

  std::cout << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}
