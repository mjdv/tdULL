#pragma once
#include <map>
#include <memory>
#include <vector>

bool IsAscending(const std::vector<int> &word) {
  for (int i = 1; i < word.size(); ++i)
    if (word[i - 1] >= word[i])
      return false;
  return true;
}

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
    if (result == children_.end())
      return nullptr;
    return &result->second;
  }
};

class SetTrie {
public:
  SetTrie() {}

  void Insert(const std::vector<int> &word) {
    assert(IsAscending(word));
    Node *node = &root_;
    for (auto n : word)
      node = node->FindOrCreateChild(n);
    node->flag_last_ = true;
  }

  Node *Search(const std::vector<int> &word) {
    assert(IsAscending(word));
    Node *node = &root_;
    for (auto n : word) {
      node = node->FindChild(n);
      if (node == nullptr)
        return nullptr;
    }

    if (node->flag_last_)
      return node;
    else
      return nullptr;
  }

  bool HasSubset(const std::vector<int> &word) {
    assert(IsAscending(word));
    return HasSubset(&root_, word, 0);
  }

  bool HasSuperset(const std::vector<int> &word) {
    assert(IsAscending(word));
    return HasSuperset(&root_, word, 0);
  }

protected:
  bool HasSubset(Node *node, const std::vector<int> &word, int idx) {
    if (node->flag_last_)
      return true;
    if (idx >= word.size())
      return false;
    bool found = false;
    Node *next_node = node->FindChild(word[idx]);
    if (next_node)
      found = HasSubset(next_node, word, idx + 1);
    if (!found)
      return HasSubset(node, word, idx + 1);
    else
      return true;
  }

  bool HasSuperset(Node *node, const std::vector<int> &word, int idx) {
    if (idx >= word.size())
      return true;
    bool found = false;
    for (auto &[num, child] : node->children_) {
      // NOTE: If we assume the children are sorted we could break.
      if (num > word[idx])
        continue;

      if (num == word[idx])
        found = HasSuperset(&child, word, idx + 1);
      else
        found = HasSuperset(&child, word, idx);

      if (found)
        break;
    }
    return found;
  }

  Node root_;
};
