#include <fstream>
#include <iostream>
#include <cassert>

#include "graph.hpp"
#include "set_trie.hpp"

#include <ctime>

SetTrie cache;

std::pair<int, int> treedepth(const SubGraph &G, int search_lbnd,
                              int search_ubnd) {

    int lower = 1, upper = G.vertices.size();

    // If you are already in the cache (exact match) we can use the lower and
    // upper from where we left off.
    auto node = cache.Search(G);
    if (node != nullptr) {
        lower = node->data.lower_bound;
        upper = node->data.upper_bound;
    }

    // If the trivial or previously found bounds suffice, we are done.
    if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
        return std::make_pair(lower, upper);
    }

    // Try to find a better lower bound by looking at cached subsets.
    for(auto subset : cache.AllSubsets(G)) {
        lower = std::max(lower, subset->data.lower_bound);
    }

    // Was this bound sufficient for our current purposes?
    if (search_ubnd <= lower || lower == upper) {
        node = cache.Insert(G);
        node->data.lower_bound = lower;
        node->data.upper_bound = upper;
        return std::make_pair(lower, upper);
    }

    // Main loop: try every vertex as root.
    // new_lower tries to find a new treedepth lower bound on this subgraph.
    int new_lower = G.vertices.size();
    for (auto v : G.vertices) {
        int search_ubnd_v = std::min(search_ubnd - 1, upper - 1);
        int search_lbnd_v = std::max(search_lbnd - 1, 1);

        int upper_v = 0;
        int lower_v = lower - 1;

        bool early_break = false;

        for (auto H : G.WithoutVertex(v).ConnectedComponents()) {
            auto pr = treedepth(H, search_lbnd_v, search_ubnd_v);
            int lower_H = pr.first;
            int upper_H = pr.second;

            upper_v = std::max(upper_v, upper_H);
            lower_v = std::max(lower_v, lower_H);

            search_lbnd_v = std::max(search_lbnd_v, lower_H);

            if (lower_H >= search_ubnd_v) {
                // This component already shows that there's no reason to
                // continue trying with vertex v.
                early_break = true;
                break;
            }
        }

        new_lower = std::min(new_lower, lower_v + 1);
        // The upper bound we found for v is only meaningful if we didn't break
        // early.
        if (upper >= upper_v + 1 && !early_break) {
            upper = upper_v + 1;
        }

        if(upper <= search_lbnd || lower == upper) {
            // Choosing root v already gives us a treedepth decomposition which
            // is good enough (a sister branch is at least this long) or
            // optimal (lower == upper) so we can use v as our root for now.
            node = cache.Insert(G);
            node->data.lower_bound = lower;
            node->data.upper_bound = upper;
            return std::make_pair(lower, upper);
        }
    }

    lower = std::max(lower, new_lower);
    node = cache.Insert(G);
    node->data.lower_bound = lower;
    node->data.upper_bound = upper;

    return std::make_pair(lower, upper);
}

std::pair<int, int> treedepth_no_subsets(const SubGraph &G, int search_lbnd,
                              int search_ubnd) {

    int lower = 1, upper = G.vertices.size();

    // If you are already in the cache (exact match) we can use the lower and
    // upper from where we left off.
    auto node = cache.Search(G);
    if (node != nullptr) {
        lower = node->data.lower_bound;
        upper = node->data.upper_bound;
    }

    // If the trivial or previously found bounds suffice, we are done.
    if (search_ubnd <= lower || search_lbnd >= upper || lower == upper) {
        return std::make_pair(lower, upper);
    }
    /*
    // Try to find a better lower bound by looking at cached subsets.
    for(auto subset : cache.AllSubsets(G)) {
        lower = std::max(lower, subset->data.lower_bound);
    }

    // Was this bound sufficient for our current purposes?
    if (search_ubnd <= lower || lower == upper) {
        node = cache.Insert(G);
        node->data.lower_bound = lower;
        node->data.upper_bound = upper;
        return std::make_pair(lower, upper);
    }*/

    // Main loop: try every vertex as root.
    // new_lower tries to find a new treedepth lower bound on this subgraph.
    int new_lower = G.vertices.size();
    for (auto v : G.vertices) {
        int search_ubnd_v = std::min(search_ubnd - 1, upper - 1);
        int search_lbnd_v = std::max(search_lbnd - 1, 1);

        int upper_v = 0;
        int lower_v = lower - 1;

        bool early_break = false;

        for (auto H : G.WithoutVertex(v).ConnectedComponents()) {
            auto pr = treedepth(H, search_lbnd_v, search_ubnd_v);
            int lower_H = pr.first;
            int upper_H = pr.second;

            upper_v = std::max(upper_v, upper_H);
            lower_v = std::max(lower_v, lower_H);

            search_lbnd_v = std::max(search_lbnd_v, lower_H);

            if (lower_H >= search_ubnd_v) {
                // This component already shows that there's no reason to
                // continue trying with vertex v.
                early_break = true;
                break;
            }
        }

        new_lower = std::min(new_lower, lower_v + 1);
        // The upper bound we found for v is only meaningful if we didn't break
        // early.
        if (upper >= upper_v + 1 && !early_break) {
            upper = upper_v + 1;
        }

        if(upper <= search_lbnd || lower == upper) {
            // Choosing root v already gives us a treedepth decomposition which
            // is good enough (a sister branch is at least this long) or
            // optimal (lower == upper) so we can use v as our root for now.
            node = cache.Insert(G);
            node->data.lower_bound = lower;
            node->data.upper_bound = upper;
            return std::make_pair(lower, upper);
        }
    }

    lower = std::max(lower, new_lower);
    node = cache.Insert(G);
    node->data.lower_bound = lower;
    node->data.upper_bound = upper;

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
        //std::cerr << "Expecting 2 arguments." << std::endl;
        return 1;
    }

    std::ifstream input;
    input.open(argv[1], std::ios::in);
    LoadGraph(input);
    input.close();

    time_t start, end;
    std::cout << "treedepth gives result:" << std::endl;

    time(&start);
    std::cout << treedepth(full_graph_as_sub, 1,
                           full_graph_as_sub.vertices.size()).second
              << std::endl;
    time(&end);
    std::cout << "Elapsed time is " << difftime(end, start) << " seconds.\n";

    cache = SetTrie();

    std::cout << "treedepth (without subset caching) gives result:" << std::endl;

    time(&start);
    std::cout << treedepth_no_subsets(full_graph_as_sub, 1,
                           full_graph_as_sub.vertices.size()).second
              << std::endl;
    time(&end);
    std::cout << "Elapsed time is " << difftime(end, start) << " seconds.\n";
}
