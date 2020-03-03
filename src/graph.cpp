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

  for (int e = 0; e < M; e++) {
    int a, b;
    stream >> a >> b;

    // Make them zero indexed.
    a--, b--;

    vertices[a].adj.emplace_back(&vertices[b]);
    vertices[b].adj.emplace_back(&vertices[a]);
  }
}

SubGraph::SubGraph(std::vector<Vertex *> &&v, std::vector<bool> &&m)
    : vertices(std::move(v)), mask(std::move(m)) {
  assert(vertices.size());
  assert(mask.size() == full_graph.N);
}
SubGraph::SubGraph() : mask(full_graph.N, false) {}

std::vector<Vertex *> SubGraph::Adj(Vertex *v) const {
  assert(mask[v->n]);
  std::vector<Vertex *> result;
  for (Vertex *child : v->adj)
    if (mask[child->n]) result.emplace_back(child);
  return result;
}

SubGraph SubGraph::WithoutVertex(Vertex *v) const {
  SubGraph result;
  result.vertices.reserve(vertices.size());
  for (Vertex *vertex : vertices)
    if (vertex != v) {
      result.vertices.emplace_back(vertex);
      result.mask[vertex->n] = true;
    }
  // Verify that we indeed remove something.
  assert(result.vertices.size() < vertices.size());
  return result;
}

std::vector<SubGraph> SubGraph::ConnectedComponents() const {
  std::vector<SubGraph> cc;
  std::vector<Vertex *> stack;

  // Initiate a DFS from all of the vertices inside this subgraph.
  for (auto root : vertices)
    if (!root->visited) {
      SubGraph component;
      component.vertices.reserve(vertices.size());

      stack.emplace_back(root);
      root->visited = true;
      while (!stack.empty()) {
        Vertex *v = stack.back();
        stack.pop_back();

        // Add this vertex to the component.
        component.vertices.emplace_back(v);
        component.mask[v->n] = true;

        // Recurse.
        for (Vertex *child : v->adj)
          if (mask[child->n] && !child->visited) {
            stack.emplace_back(child);
            child->visited = true;
          }
      }
      cc.emplace_back(std::move(component));
    }

  // Reset the visited field.
  for (auto vtx : vertices) vtx->visited = false;
  return cc;
}

void LoadGraph(std::istream &stream) {
  full_graph = Graph(stream);

  // Create the subgraph.
  std::vector<bool> mask(full_graph.N, true);
  std::vector<Vertex *> vertices;
  vertices.reserve(full_graph.N);
  for (Vertex &vertex : full_graph.vertices) vertices.emplace_back(&vertex);
  full_graph_as_sub = SubGraph(std::move(vertices), std::move(mask));

  std::cout << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}
