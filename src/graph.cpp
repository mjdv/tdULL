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

    adj[a].emplace_back(&vertices[b]);
    adj[b].emplace_back(&vertices[a]);
  }
}

SubGraph::SubGraph() : mask(full_graph.N, false) {}

const std::vector<Vertex *> &SubGraph::Adj(Vertex *v) const {
  assert(mask[v->n]);
  return adj.at(v);
}

SubGraph SubGraph::WithoutVertex(Vertex *v) const {
  assert(mask[v->n]);
  SubGraph result;
  result.mask = mask;
  result.mask[v->n] = false;

  result.vertices.reserve(vertices.size());
  for (Vertex *vertex : vertices) {
    if (vertex == v) continue;
    result.vertices.emplace_back(vertex);

    // Find the (trimmed) adjacency list for vertex.
    std::vector<Vertex *> nghbrs;
    nghbrs.reserve(Adj(vertex).size());
    for (Vertex *nghb : Adj(vertex))
      if (nghb != v) nghbrs.emplace_back(nghb);

    // Insert this updated adjacency list.
    result.M += nghbrs.size();
    result.adj.emplace(vertex, std::move(nghbrs));
  }

  assert(result.vertices.size() == result.adj.size());

  // Verify that we indeed remove something.
  assert(result.vertices.size() < vertices.size());
  return result;
}

std::vector<SubGraph> SubGraph::ConnectedComponents() const {
  std::vector<SubGraph> cc;
  static std::vector<Vertex *> stack;

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = vertices.size();
  for (auto root : vertices) {
    if (vertices_left == 0) break;
    if (!root->visited) {
      SubGraph component;
      component.vertices.reserve(vertices_left);

      stack.emplace_back(root);
      root->visited = true;
      while (!stack.empty()) {
        Vertex *v = stack.back();
        stack.pop_back();

        // Add this vertex and its adjacency to the component.
        component.vertices.emplace_back(v);
        component.mask[v->n] = true;
        component.adj.emplace(v, Adj(v));
        component.M += Adj(v).size();
        vertices_left--;

        // Recurse.
        for (Vertex *nghb : Adj(v))
          if (!nghb->visited) {
            stack.emplace_back(nghb);
            nghb->visited = true;
          }
      }
      cc.emplace_back(std::move(component));
    }
  }

  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return cc;
}

void LoadGraph(std::istream &stream) {
  full_graph = Graph(stream);

  // Create the subgraph.
  full_graph_as_sub.mask = std::vector<bool>(full_graph.N, true);
  for (Vertex &vertex : full_graph.vertices) {
    full_graph_as_sub.vertices.emplace_back(&vertex);
    full_graph_as_sub.adj.emplace(&vertex, full_graph.adj.at(vertex.n));
  }
  full_graph_as_sub.M = full_graph.M;

  std::cout << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}
