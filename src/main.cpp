#include "graph.hpp"
#include "set_trie.hpp"
#include <iostream>
#include <fstream>

//SetTrie cache;

std::pair<int, int> treedepth(SubGraph G, int search_lbnd, int search_ubnd) {
    if(search_ubnd <= 1 || search_lbnd >= G.vertices_.size()) {
        //cache.Insert(...) skip caching for now
        return std::make_pair(1, G.vertices_.size());
    }

    
}

int treedepth_trivial(SubGraph G) {
    if(G.vertices_.size() == 1)
        return 1;
    
    int td = G.vertices_.size();
    for(auto v : G.vertices_) {
        int td_for_this_root = 0;
        
        for(auto cc : G.WithoutVertex(v).ConnectedComponents()) {
            td_for_this_root = std::max(td_for_this_root, treedepth_trivial(cc) + 1);
        }
        td = std::min(td, td_for_this_root);
    }
    return td;
}

int main(int argc, char **argv) {
    if(argc != 3) {
        std::cerr << "Expecting 2 arguments." << std::endl;
        return 1;
    }

    std::ifstream input;
    input.open(argv[1], std::ios::in);
    LoadGraph(input);
    input.close();

    std::cout << treedepth_trivial(full_graph_as_sub) << std::endl;
}
