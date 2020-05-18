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
  virtual bool HasNext() const = 0;
  virtual std::vector<Separator> Next(int k = 10000) = 0;
  virtual ~SeparatorGenerator() = default;
};

class SeparatorGeneratorDirected : public SeparatorGenerator {
 public:
  SeparatorGeneratorDirected(const Graph &G, const int source);

  bool HasNext() const override { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000) override;

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // Reference to the graph for which we are generating separators.
  const Graph &G;
  int source;

  // In done we keep the seperators we have already enqueued, to make sure
  // they aren't processed again. In queue we keep all the ones we have
  // generated, but which we have not yet used to generate new ones.
  std::queue<std::pair<std::vector<int>, std::vector<bool>>> queue;
  std::unordered_set<std::vector<bool>> done;

  // In buffer we will keep all the (fully minimal) generated separators.
  std::vector<Separator> buffer;

  // Shared datatypes.
  std::vector<bool> in_nbh;
  std::vector<bool> sep_mask;
};

class SeparatorGeneratorUndirected : public SeparatorGenerator {
 public:
  SeparatorGeneratorUndirected(const Graph &G);

  bool HasNext() const override { return !queue.empty(); }
  std::vector<Separator> Next(int k = 10000) override;

  void clear() {
    done.clear();
    queue = {};
  }

 protected:
  // Reference to the graph for which we are generating separators.
  const Graph &G;

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
};
