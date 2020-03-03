#include <fstream>
#include <iostream>

#include "graph.hpp"
#include "set_trie.hpp"

SetTrie cache;

void CacheInsert(const SubGraph &G, int lower, int upper, int root) {
    auto node = cache.Insert(G);
    node->data.lower_bound = lower;
    node->data.upper_bound = upper;
    node->data.root = root;
}

std::pair<int, int> treedepth(const SubGraph &G, int search_lbnd,
                              int search_ubnd) {
    /*if (search_ubnd <= 1 || search_lbnd >= G.vertices.size()) {
        CacheInsert(G, 1, G.vertices.size(), G.vertices[0]->n);
        return std::make_pair(1, G.vertices.size());
    }*/

    int lower = 1, upper = G.vertices.size();
    int current_root = G.vertices[0]->n;

    auto node = cache.Search(G);
    if(node != nullptr) {
        lower = node->data.lower_bound;
        upper = node->data.upper_bound;
        current_root = node->data.root;
    }

    if (search_ubnd <= lower || search_lbnd >= upper) {
        CacheInsert(G, lower, upper, current_root);
        return std::make_pair(lower, upper);
    }

    int new_lower = 0;
    auto sorted_vertices = G.vertices;
    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
    [&](Vertex *v1, Vertex *v2) {
        return G.Adj(v1).size() > G.Adj(v2).size();
    });
    assert(G.Adj(sorted_vertices[0]).size() >= G.Adj(sorted_vertices[1]).size());
    for(auto v : G.vertices) {
        int search_ubnd_v = std::min(search_ubnd - 1, upper - 1);
        int search_lbnd_v = std::max(search_lbnd - 1, lower - 1);

        int upper_from_v = 0;
        int lower_from_v = lower - 1;

        for (auto H : G.WithoutVertex(v).ConnectedComponents()) {
            auto pr = treedepth(H, search_lbnd - 1, search_ubnd_v);
            int lower_H = pr.first;
            int upper_H = pr.second;

            if (lower_H >= search_ubnd_v) break;

            upper_from_v = std::max(upper_from_v, upper_H);
            lower_from_v = std::min(lower_from_v, lower_H);

            search_lbnd_v = std::max(search_lbnd_v, lower_H);
        }

        new_lower = std::max(new_lower, lower_from_v + 1);
        if (upper <= upper_from_v) {
            upper = upper_from_v;
            current_root = v->n;
        }

        if(search_lbnd >= upper || lower == upper) {
            CacheInsert(G, lower, upper, current_root);
            return std::make_pair(lower, upper);
        }
    }

    lower = std::min(lower, new_lower);

    CacheInsert(G, lower, upper, current_root);
    return std::make_pair(lower, upper);
}

int treedepth_trivial(const SubGraph &G) {
    if (G.vertices.size() == 1) return 1;

    auto node = cache.Search(G);
    if (node != nullptr) {
        return node->data.upper_bound;
    }

    int td = G.vertices.size();
    for (auto v : G.vertices) {
        int td_for_this_root = 0;

        for (auto cc : G.WithoutVertex(v).ConnectedComponents()) {
            td_for_this_root = std::max(td_for_this_root, treedepth_trivial(cc) + 1);
        }
        td = std::min(td, td_for_this_root);
    }

    node = cache.Insert(G);
    node->data.upper_bound = td;
    return td;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Expecting 2 arguments." << std::endl;
        return 1;
    }

    std::ifstream input;
    input.open(argv[1], std::ios::in);
    LoadGraph(input);
    input.close();

    std::cout << treedepth(full_graph_as_sub, 1, full_graph_as_sub.vertices.size()).second << std::endl;
}
