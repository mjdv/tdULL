#include "set_trie.hpp"

#include <iostream>

bool IsAscending(const std::vector<int> &word) {
  for (int i = 1; i < word.size(); ++i)
    if (word[i - 1] >= word[i]) return false;
  return true;
}

void SetTrie::Insert(const std::vector<int> &word) {
  assert(IsAscending(word));
  Node *node = &root_;
  for (auto n : word) node = node->FindOrCreateChild(n);
  node->flag_last_ = true;
}

Node *SetTrie::Search(const std::vector<int> &word) {
  assert(IsAscending(word));
  Node *node = &root_;
  for (auto n : word) {
    node = node->FindChild(n);
    if (node == nullptr) return nullptr;
  }

  if (node->flag_last_)
    return node;
  else
    return nullptr;
}

bool SetTrie::HasSubset(const std::vector<int> &word) {
  assert(IsAscending(word));
  std::stack<std::pair<Node *, int>, std::vector<std::pair<Node *, int>>> stack;
  stack.emplace(&root_, 0);

  while (!stack.empty()) {
    auto [node, idx] = stack.top();
    stack.pop();
    if (node->flag_last_) return true;
    if (idx >= word.size()) continue;
    stack.emplace(node, idx + 1);
    Node *next_node = node->FindChild(word[idx]);
    if (next_node) stack.emplace(next_node, idx + 1);
  }
  return false;
}

bool SetTrie::HasSuperset(const std::vector<int> &word) {
  assert(IsAscending(word));
  std::stack<std::pair<Node *, int>, std::vector<std::pair<Node *, int>>> stack;
  stack.emplace(&root_, 0);

  while (!stack.empty()) {
    auto [node, idx] = stack.top();
    stack.pop();
    if (idx >= word.size()) return true;
    for (auto &[num, child] : node->children_) {
      // NOTE: If we assume the children are sorted we could break.
      if (num > word[idx]) continue;

      if (num == word[idx])
        stack.emplace(&child, idx + 1);
      else
        stack.emplace(&child, idx);
    }
  }
  return false;
}
