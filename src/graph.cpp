#include "graph.hpp"

Graph full_graph;
SubGraph full_graph_as_sub;

Graph::Graph(std::istream &stream) {
  std::string str;
  stream >> str;
  assert(str == "p");
  stream >> str;
  assert(str == "tdp");
  stream >> N_ >> M_;

  // Create the vector of vertices.
  vertices_.reserve(N_);
  for (int v = 0; v < N_; v++) vertices_.emplace_back(v);

  for (int e = 0; e < M_; e++) {
    int a, b;
    stream >> a >> b;

    // Make them zero indexed.
    a--, b--;

    vertices_[a].adj_.emplace_back(&vertices_[b]);
    vertices_[b].adj_.emplace_back(&vertices_[a]);
  }
}

SubGraph::SubGraph(std::vector<Vertex *> &&vertices, std::vector<bool> &&mask)
    : vertices_(std::move(vertices)), mask_(std::move(mask)) {
  assert(vertices_.size());
  assert(mask_.size() == full_graph.N_);
}
SubGraph::SubGraph() : mask_(full_graph.N_, false) {}

SubGraph SubGraph::WithoutVertex(Vertex *v) {
  SubGraph result;
  result.vertices_.reserve(vertices_.size());
  for (Vertex *vertex : vertices_)
    if (vertex != v) {
      result.vertices_.emplace_back(vertex);
      result.mask_[vertex->n_] = true;
    }
  // Verify that we indeed remove something.
  assert(result.vertices_.size() < vertices_.size());
  return result;
}

std::vector<SubGraph> SubGraph::ConnectedComponents() {
  std::vector<SubGraph> cc;
  std::vector<Vertex *> stack;

  // Initiate a DFS from all of the vertices inside this subgraph.
  for (auto root : vertices_)
    if (!root->visited_) {
      SubGraph component;
      component.vertices_.reserve(vertices_.size());

      stack.emplace_back(root);
      root->visited_ = true;
      while (!stack.empty()) {
        Vertex *v = stack.back();
        stack.pop_back();

        // Add this vertex to the component.
        component.vertices_.emplace_back(v);
        component.mask_[v->n_] = true;

        // Recurse.
        for (Vertex *child : v->adj_)
          if (mask_[child->n_] && !child->visited_) {
            stack.emplace_back(child);
            child->visited_ = true;
          }
      }
      cc.emplace_back(std::move(component));
    }

  // Reset the visited field.
  for (auto vtx : vertices_) vtx->visited_ = false;
  return cc;
}

void LoadGraph(std::istream &stream) {
  full_graph = Graph(stream);

  // Create the subgraph.
  std::vector<bool> mask(full_graph.N_, true);
  std::vector<Vertex *> vertices;
  vertices.reserve(full_graph.N_);
  for (Vertex &vertex : full_graph.vertices_) vertices.emplace_back(&vertex);
  full_graph_as_sub = SubGraph(std::move(vertices), std::move(mask));

  std::cout << "Initalized a graph having " << full_graph.N_
            << " vertices with " << full_graph.M_ << " edges. " << std::endl;
}
