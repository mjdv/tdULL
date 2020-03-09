#pragma once
#include <cassert>
#include <climits>
#include <map>
#include <memory>
#include <stack>
#include <vector>

struct Node {
  // This is the data that will be stored inside the Set Trie.
  int16_t upper_bound = INT16_MAX;
  int16_t lower_bound = 0;
  int16_t root = -1;
  bool flag_last = false;

  // TODO: Is std::map the best datastructure? Sorting does help.
  std::map<int16_t, Node> children;

  // Returns pointer to child, and nullptr if it doesn't exist.
  Node *FindChild(int16_t n) {
    auto result = children.find(n);
    if (result == children.end()) return nullptr;
    return &result->second;
  }

  // Find a child with the given letter, creates one if it doesn't yet exist.
  Node *FindOrCreateChild(int16_t n) {
    return &children.emplace(n, Node()).first->second;
  }
};

class SetTrie {
 public:
  std::pair<Node *, bool> Insert(const std::vector<int16_t> &word);
  Node *Search(const std::vector<int16_t> &word);

  bool HasSubset(const std::vector<int16_t> &word);
  bool HasSuperset(const std::vector<int16_t> &word);

  std::vector<Node *> AllSubsets(const std::vector<int16_t> &word);
  std::vector<Node *> AllSupersets(const std::vector<int16_t> &word);

 protected:
  Node root_;
};
