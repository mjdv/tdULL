#pragma once
#include <cassert>

#include "graph.hpp"

struct Separator {
  std::vector<int> vertices;
  std::pair<int, int> largest_component;
  bool fully_minimal = false;

  Separator(const Graph &G, const std::vector<int> &vertices);
};

class SeparatorGenerator {
 public:
  virtual bool HasNext() = 0;
  virtual std::vector<Separator> Next(int k = 10000) = 0;
  virtual ~SeparatorGenerator() = default;
  SeparatorGenerator(const Graph &G);

 protected:
  // Reference to the graph for which we are generating separators.
  const Graph &G;

  // Stores the separators that we have already visited.
  std::unordered_set<std::vector<bool>> done;

  // In buffer we will keep all the (fully minimal) generated separators.
  std::vector<Separator> buffer;

  // Shared datatypes.
  std::vector<bool> in_nbh;
  std::vector<bool> sep_mask;
};


class SeparatorGeneratorDirected : public SeparatorGenerator {
 public:
  SeparatorGeneratorDirected(const Graph &G, const int source);

  bool HasNext() override { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000) override;

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // The source from which we generate separators them.
  int source;

  // In done we keep the seperators we have already enqueued, to make sure
  // they aren't processed again. In queue we keep all the ones we have
  // generated, but which we have not yet used to generate new ones.
  // Additionally the queue contains a bitmask encoding the component in the
  // complement of source: this saves us doing a DFS every time.
  std::queue<std::pair<std::vector<int>, std::vector<bool>>> queue;
};

class SeparatorGeneratorUndirected : public SeparatorGenerator {
 public:
  SeparatorGeneratorUndirected(const Graph &G);

  bool HasNext() override { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000) override;

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // In done we keep the seperators we have already enqueued, to make sure
  // they aren't processed again. In queue we keep all the ones we have
  // generated, but which we have not yet used to generate new ones.
  std::queue<std::vector<int>> queue;
};

class SeparatorGenerator2Core : public SeparatorGenerator {
 public:
  SeparatorGenerator2Core(const Graph &G);

  bool HasNext() override;

  std::vector<Separator> Next(int k = 10000) override;

  void clear() {
    done.clear();
    two_core_gen.clear();
  }

 protected:
  std::vector<int> G_to_two_core;
  std::vector<int> two_core_to_G;
  std::vector<bool> two_core_adj;
  const Graph two_core;
  SeparatorGeneratorUndirected two_core_gen;
};
