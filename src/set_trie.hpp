#pragma once
#include <cassert>
#include <map>
#include <memory>
#include <stack>
#include <vector>

struct Node {
  bool flag_last_ = false;

  // TODO: Is std::map the best datastructure?
  std::map<int, Node> children_;

  // Find a child with the given letter, creates one if it doesn't yet exist.
  Node *FindOrCreateChild(int n) {
    return &children_.emplace(n, Node()).first->second;
  }

  // Returns pointer to child, and nullptr if it doesn't exist.
  Node *FindChild(int n) {
    auto result = children_.find(n);
    if (result == children_.end()) return nullptr;
    return &result->second;
  }
};

class SetTrie {
 public:
  SetTrie() {}

  void Insert(const std::vector<int> &word);
  Node *Search(const std::vector<int> &word);

  bool HasSubset(const std::vector<int> &word);
  bool HasSuperset(const std::vector<int> &word);

 protected:
  Node root_;
};
