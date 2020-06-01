#pragma once
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <parallel_hashmap/phmap.h>

#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <unordered_set>

#include "graph.hpp"

struct Separator {
  std::vector<int> vertices;
  std::pair<int, int> largest_component;

  // NOTE: This is either a fully minimal separator, or it is not fully minimal
  // but the non-minimality comes from leaves being cut off.
  bool fully_minimal = false;
  size_t num_components;

  Separator(const Graph &G, const std::vector<int> &vertices);
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

  // Reference to the graph for which we are generating separators.
  const Graph &G;

  // In done we keep the seperators we have already enqueued, to make sure
  // they aren't processed again. In queue we keep all the ones we have
  // generated, but which we have not yet used to generate new ones.
  std::queue<std::vector<int>> queue;
  struct BitsetHash {
    template <typename Block, typename Allocator>
    inline size_t operator()(
        const boost::dynamic_bitset<Block, Allocator> &a) const BOOST_NOEXCEPT {
      std::size_t res = boost::hash_value(a.m_num_bits);
      boost::hash_combine(res, a.m_bits);
      return res;
    }
  };
  phmap::flat_hash_set<boost::dynamic_bitset<>, BitsetHash> done;

  // In buffer we will keep all the (fully minimal) generated separators.
  std::vector<Separator> buffer;

  // Shared datatypes.
  std::vector<bool> in_nbh;
  boost::dynamic_bitset<> sep_mask;
};
