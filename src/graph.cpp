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

    adj[a].emplace_back(&vertices[b]);
    adj[b].emplace_back(&vertices[a]);
  }
}

SubGraph::SubGraph() : mask(full_graph.N, false) {}

const std::vector<Vertex *> &SubGraph::Adj(Vertex *v) const {
  assert(mask[v->n]);
  return adj.at(v);
}

std::vector<SubGraph> SubGraph::WithoutVertex(Vertex *w) const {
  assert(mask[w->n]);
  std::vector<SubGraph> cc;
  static std::vector<Vertex *> stack;

  // Initiate a DFS from all of the vertices inside this subgraph.
  int vertices_left = vertices.size();
  for (auto root : vertices) {
    if (vertices_left == 0) break;
    if (root == w) continue;
    if (!root->visited) {
      SubGraph component;
      component.vertices.reserve(vertices_left);

      stack.emplace_back(root);
      root->visited = true;
      while (!stack.empty()) {
        Vertex *v = stack.back();
        stack.pop_back();

        // Insert this vertex into the component.
        component.vertices.emplace_back(v);
        component.mask[v->n] = true;
        vertices_left--;

        // Find the (trimmed) adjacency list for vertex.
        std::vector<Vertex *> nghbrs;
        nghbrs.reserve(Adj(v).size());
        for (Vertex *nghb : Adj(v))
          if (nghb != w) nghbrs.emplace_back(nghb);

        // Insert this adjacency list into the component.
        component.M += nghbrs.size();
        component.max_degree = std::max(component.max_degree, nghbrs.size());
        auto [it, inserted] = component.adj.emplace(v, std::move(nghbrs));
        assert(inserted);

        // Recurse.
        for (Vertex *nghb : it->second)
          if (!nghb->visited) {
            stack.emplace_back(nghb);
            nghb->visited = true;
          }
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

std::vector<Vertex *> SubGraph::Bfs(Vertex *root) const {
  assert(mask[root->n]);

  std::vector<Vertex *> result;
  result.reserve(vertices.size());

  static std::queue<Vertex *> queue;
  queue.push(root);
  root->visited = true;
  while (!queue.empty()) {
    Vertex *v = queue.front();
    queue.pop();
    result.emplace_back(v);
    for (Vertex *nghb : Adj(v))
      if (!nghb->visited) {
        queue.push(nghb);
        nghb->visited = true;
      }
  }
  // Reset the visited field.
  for (auto v : result) v->visited = false;
  return result;
}

void LoadGraph(std::istream &stream) {
  full_graph = Graph(stream);

  // Create the subgraph.
  full_graph_as_sub.M = full_graph.M;
  full_graph_as_sub.mask = std::vector<bool>(full_graph.N, true);
  for (Vertex &vertex : full_graph.vertices) {
    full_graph_as_sub.vertices.emplace_back(&vertex);
    full_graph_as_sub.adj.emplace(&vertex, full_graph.adj.at(vertex.n));
    full_graph_as_sub.max_degree = std::max(full_graph_as_sub.max_degree,
                                            full_graph.adj.at(vertex.n).size());
  }

  std::cout << "Initalized a graph having " << full_graph.N << " vertices with "
            << full_graph.M << " edges. " << std::endl;
}
