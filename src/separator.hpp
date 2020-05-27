#pragma once
#include <memory>

#include "graph.hpp"
#include "nauty.hpp"

struct Separator {
  std::vector<int> vertices;
  std::pair<int, int> largest_component;

  // NOTE: This is either a fully minimal separator, or it is not fully minimal
  // but the non-minimality comes from leaves being cut off.
  bool fully_minimal = false;

  Separator(const Graph &G, std::vector<int> &&vertices);
};

class SeparatorGenerator {
 public:
  SeparatorGenerator(const Graph &G);

  bool HasNext() const { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000);

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // Reference to the graph for which we are generating separators.
  const Graph &G_orig;

  // Contracted graph.
  Graph G;
  std::vector<std::vector<int>> vertices_original;

  // In done we keep the seperators we have already enqueued, to make sure
  // they aren't processed again. In queue we keep all the ones we have
  // generated, but which we have not yet used to generate new ones.
  std::queue<std::vector<int>> queue;
  std::unordered_set<std::vector<bool>> done;

  // In buffer we will keep all the (fully minimal) generated separators.
  std::vector<Separator> buffer;

  // Shared datatypes.
  std::vector<bool> in_nbh;
  std::vector<bool> sep_mask;
  std::vector<std::vector<int>> in_automorphisms;
  std::unique_ptr<Nauty> nauty;
};
