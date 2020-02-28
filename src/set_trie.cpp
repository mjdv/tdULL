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

// Recursive impl.
bool HasSubset(Node *node, const std::vector<int> &word, int idx) {
  if (node->flag_last_) return true;
  if (idx >= word.size()) return false;
  bool found = false;
  Node *next_node = node->FindChild(word[idx]);
  if (next_node) found = HasSubset(next_node, word, idx + 1);
  if (!found)
    return HasSubset(node, word, idx + 1);
  else
    return true;
}

bool SetTrie::HasSubset(const std::vector<int> &word) {
  assert(IsAscending(word));
  return ::HasSubset(&root_, word, 0);
}

// Recursive impl.
bool HasSuperset(Node *node, const std::vector<int> &word, int idx) {
  if (idx >= word.size()) return true;
  bool found = false;
  for (auto &[num, child] : node->children_) {
    // NOTE: If we assume the children are sorted we could break.
    if (num > word[idx]) continue;

    if (num == word[idx])
      found = HasSuperset(&child, word, idx + 1);
    else
      found = HasSuperset(&child, word, idx);

    if (found) break;
  }
  return found;
}

bool SetTrie::HasSuperset(const std::vector<int> &word) {
  assert(IsAscending(word));
  return ::HasSuperset(&root_, word, 0);
}
