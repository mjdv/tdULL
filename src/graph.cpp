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
  assert(result.M == total_edges/2);
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
  assert(result.M == total_edges/2);
  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return result;
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
