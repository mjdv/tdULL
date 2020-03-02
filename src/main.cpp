#include "graph.hpp"
#include "set_trie.hpp"
#include <iostream>

//SetTrie cache;

std::pair<int, int> treedepth(SubGraph G, int search_lbnd, int search_ubnd) {
    if(search_ubnd <= 1 || search_lbnd >= G.vertices_.size()) {
        //cache.Insert(...) skip caching for now
        return std::make_pair(1, G.vertices_.size());
    }
}

int main() {
    LoadGraph(std::cin);
}
