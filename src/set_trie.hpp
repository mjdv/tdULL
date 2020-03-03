#pragma once
#include <cassert>
#include <map>
#include <memory>
#include <stack>
#include <vector>
#include <climits>

// The data that we will store inside the Set Trie.
struct Data {
  int upper_bound = INT_MAX;
  int lower_bound = 0;
  int root = -1;
};

struct Node {
  Data data;  // Store the data and make it publicly available.
  int n;      // The number of this node.

  // TODO: Is std::map the best datastructure? Sorting does help.
  std::map<int, Node> children;
  Node *parent;
  bool flag_last = false;

  Node(int n, Node *parent, const Data &data)
      : data(data), n(n), parent(parent) {}
  Node(int n, Node *parent) : Node(n, parent, {}) {}

  // Returns the word associated to this node.
  std::vector<int> Word() const {
    std::vector<int> word;
    auto node = this;
    while (node && node->parent) {
      word.push_back(node->n);
      node = node->parent;
    }
    return {word.rbegin(), word.rend()};
  }

  // Returns pointer to child, and nullptr if it doesn't exist.
  Node *FindChild(int n) {
    auto result = children.find(n);
    if (result == children.end()) return nullptr;
    return &result->second;
  }

  // Find a child with the given letter, creates one if it doesn't yet exist.
  Node *FindOrCreateChild(int n) {
    return &children.emplace(n, Node{n, this}).first->second;
  }
};

class SetTrie {
 public:
  SetTrie() : root_(-1, nullptr) {}

  Node *Insert(const std::vector<int> &word);
  Node *Search(const std::vector<int> &word);

  bool HasSubset(const std::vector<int> &word);
  bool HasSuperset(const std::vector<int> &word);

  std::vector<Node *> AllSubsets(const std::vector<int> &word);
  std::vector<Node *> AllSupersets(const std::vector<int> &word);

 protected:
  Node root_;
};
