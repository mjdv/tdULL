#pragma once
#include <cassert>
#include <map>
#include <memory>
#include <stack>
#include <vector>

struct Node {
  bool flag_last_ = false;
  int n_;  // The letter of this node.

  // TODO: Is std::map the best datastructure? Sorting does help.
  std::map<int, Node> children_;
  Node *parent_;

  Node(int n, Node *parent) : n_(n), parent_(parent) {}

  // Find a child with the given letter, creates one if it doesn't yet exist.
  Node *FindOrCreateChild(int n) {
    return &children_.emplace(n, Node{n, this}).first->second;
  }

  // Returns pointer to child, and nullptr if it doesn't exist.
  Node *FindChild(int n) {
    auto result = children_.find(n);
    if (result == children_.end()) return nullptr;
    return &result->second;
  }

  // Returns the word associated to this node.
  std::vector<int> Word() const {
    std::vector<int> word;
    auto node = this;
    while (node) {
      word.push_back(node->n_);
      node = node->parent_;
    }
    return {word.rbegin(), word.rend()};
  }
};

class SetTrie {
 public:
  SetTrie() : root_(-1, nullptr) {}

  void Insert(const std::vector<int> &word);
  Node *Search(const std::vector<int> &word);

  bool HasSubset(const std::vector<int> &word);
  bool HasSuperset(const std::vector<int> &word);

  std::vector<Node *> AllSubsets(const std::vector<int> &word);

 protected:
  Node root_;
};
