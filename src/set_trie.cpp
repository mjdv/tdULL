#include "set_trie.hpp"

#include <climits>
#include <iostream>

bool IsAscending(const std::vector<int> &word) {
  for (int i = 1; i < word.size(); ++i)
    if (word[i - 1] >= word[i]) return false;
  return true;
}

Node *SetTrie::Insert(const std::vector<int> &word) {
  assert(IsAscending(word));
  Node *node = &root_;
  for (auto n : word) node = node->FindOrCreateChild(n);
  node->flag_last = true;
  return node;
}

Node *SetTrie::Search(const std::vector<int> &word) {
  assert(IsAscending(word));
  Node *node = &root_;
  for (auto n : word) {
    node = node->FindChild(n);
    if (node == nullptr) return nullptr;
  }

  if (node->flag_last)
    return node;
  else
    return nullptr;
}

// Recursive impl.
bool HasSubset(Node *node, const std::vector<int> &word, int idx) {
  if (node->flag_last) return true;
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
  for (auto &[num, child] : node->children) {
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

//
void AllSubsets(Node *node, const std::vector<int> &word, int idx,
                std::vector<Node *> &result) {
  if (node->flag_last) result.push_back(node);
  if (idx >= word.size()) return;
  for (auto &[num, child] : node->children)
    for (int j = idx; j < word.size(); ++j) {
      // TODO: If children is sorted we can break earlier.
      // if (num < word[j]) break;

      if (num == word[j]) {
        AllSubsets(&child, word, j, result);
        break;
      }
    }
}

std::vector<Node *> SetTrie::AllSubsets(const std::vector<int> &word) {
  assert(IsAscending(word));
  std::vector<Node *> result;
  ::AllSubsets(&root_, word, 0, result);
  return result;
}

void AllSupersets(Node *node, const std::vector<int> &word, int idx,
                  std::vector<Node *> &result) {
  if (idx >= word.size()) result.push_back(node);
  int current_letter = idx < word.size() ? word[idx] : INT_MAX;
  for (auto &[num, child] : node->children) {
    // TODO: If children is sorted we can break.
    if (num > current_letter) continue;

    if (num == current_letter)
      AllSupersets(&child, word, idx + 1, result);
    else
      AllSupersets(&child, word, idx, result);
  }
}

std::vector<Node *> SetTrie::AllSupersets(const std::vector<int> &word) {
  assert(IsAscending(word));
  std::vector<Node *> result;
  ::AllSupersets(&root_, word, 0, result);
  return result;
}
