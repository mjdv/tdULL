#pragma once
#include <cassert>
#include <climits>
#include <map>
#include <memory>
#include <stack>
#include <vector>

struct Node {
  // This is the data that will be stored inside the Set Trie.
  int upper_bound = INT_MAX;
  int lower_bound = 0;
  int root = -1;

  // TODO: Is std::map the best datastructure? Sorting does help.
  std::map<int, Node> children;
  bool flag_last = false;

  // Returns pointer to child, and nullptr if it doesn't exist.
  Node *FindChild(int n) {
    auto result = children.find(n);
    if (result == children.end()) return nullptr;
    return &result->second;
  }

  // Find a child with the given letter, creates one if it doesn't yet exist.
  Node *FindOrCreateChild(int n) {
    return &children.emplace(n, Node()).first->second;
  }

  size_t sizeRecur() {
    size_t result = children.size();
    for (auto &[_, child] : children) result += child.sizeRecur();
    return result;
  }
};

class SetTrie {
 public:
  std::pair<Node *, bool> Insert(const std::vector<int> &word);
  Node *Search(const std::vector<int> &word);

  bool HasSubset(const std::vector<int> &word);
  bool HasSuperset(const std::vector<int> &word);

  std::vector<Node *> AllSubsets(const std::vector<int> &word);
  std::vector<Node *> AllSupersets(const std::vector<int> &word);

  size_t size() { return root_.sizeRecur(); }

 protected:
  Node root_;
};
