#include "set_trie.hpp"

int main() {
    SetTrie cache;
    cache.Insert({1, 20, 50});
    cache.Insert({5, 8, 15});

    assert(cache.Search({1, 20, 50}));
    assert(cache.Search({5, 8, 15}));
    assert(!cache.Search({5, 8, 20}));

    assert(cache.HasSubset({5, 8, 15, 20}));
    assert(cache.HasSubset({5, 8, 10, 15}));
    assert(!cache.HasSubset({5, 8, 10}));

    assert(cache.HasSuperset({}));
    assert(cache.HasSuperset({5}));
    assert(cache.HasSuperset({5, 15}));
    assert(!cache.HasSuperset({5, 20}));
    return 0;
}
