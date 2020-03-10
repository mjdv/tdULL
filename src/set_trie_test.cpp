#include "set_trie.hpp"

#include <assert.h>

#include <climits>
#include <iostream>

int main() {
  SetTrie cache;
  cache.Insert({1, 20, 50});
  cache.Insert({5, 8, 15});

  assert(cache.Search({1, 20, 50}));
  assert(cache.Search({5, 8, 15}));
  assert(!cache.Search({5, 8, 20}));

  assert(cache.HasSubset({5, 8, 15}));
  assert(cache.HasSubset({5, 8, 15, 20}));
  assert(cache.HasSubset({5, 8, 10, 15}));
  assert(!cache.HasSubset({5, 8, 10}));
  assert(!cache.HasSubset({}));

  assert(cache.HasSuperset({}));
  assert(cache.HasSuperset({5}));
  assert(cache.HasSuperset({5, 15}));
  assert(!cache.HasSuperset({5, 20}));

  assert(cache.BigSubsets({1, 5, 8, 15, 20, 50}, 3).size() == 2);
  assert(cache.BigSubsets({1, 5, 8, 15, 20, 50}, 4).size() == 2);
  assert(cache.BigSubsets({1, 5, 8, 15, 20, 50}, 2).size() == 0);

  cache.Insert({5});
  cache.Insert({5, 8});
  cache.Insert({8, 15});

  std::cout << "Subsets{5, 7, 8, 9, 15}:" << std::endl;
  for (auto node : cache.AllSubsets({5, 7, 8, 9, 15})) {
    std::cout << "\t{";
    // for (auto l : node->Word()) std::cout << l << " ";
    std::cout << "}" << std::endl;
  }

  for (int gap = 2; gap < 4; gap++) {
    std::cout << "Subsets{5, 7, 8, 9, 15} with max gap " << gap << ":"
              << std::endl;
    for (auto node : cache.BigSubsets({5, 7, 8, 9, 15}, gap)) {
      std::cout << "\t{";
      //for (auto l : node->Word()) std::cout << l << " ";
      std::cout << "}" << std::endl;
    }
  }

  assert(cache.HasSuperset({8}));
  assert(cache.HasSuperset({5, 8, 15}));
  assert(cache.HasSuperset({8, 15}));
  std::cout << "Supersets{8}:" << std::endl;
  for (auto node : cache.AllSupersets({8})) {
    std::cout << "\t{";
    // for (auto l : node->Word()) std::cout << l << " ";
    std::cout << "}" << std::endl;
  }

  return 0;
}
