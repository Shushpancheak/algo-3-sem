/*
 * C. Шаблон с ?
 * Все языки	make2
 * Ограничение времени	3 секунды	80 секунд
 * Ограничение памяти	12Mb	400Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Шаблон поиска задан строкой длины m, в которой кроме
 * обычных символов могут встречаться символы “?”. Найти
 * позиции всех вхождений шаблона в тексте длины n. Каждое
 * вхождение шаблона предполагает, что все обычные символы
 * совпадают с соответствующими из текста, а вместо символа
 * “?” в тексте встречается произвольный символ.
 * Время работы - O(n + m + Z), где Z - общее -число
 * вхождений подстрок шаблона “между вопросиками”
 * в исходном тексте. m ≤ 5000, n ≤ 2000000.
 *
 * T(n, m, Z) = O(n + m + Z)
 * M(n, m, Z) = O(n + m + Z)
 */
#include <vector>
#include <string>
#include <iostream>

const unsigned ALPHABET_SIZE = 26;
const char AHO_KORASICK_ROOT_CHAR = 'Z';
const char AHO_KORASICK_DEFAULT_CHAR  = '\0';

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - -DECLARATIONS - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

struct Node {
  Node();
  bool IsRoot() const;
  void MakeRoot();
  void MakeTerminal(unsigned pattern_num);

  Node* son_node[ALPHABET_SIZE];
  Node* transit [ALPHABET_SIZE];
  Node* parent;
  Node* suff_link;
  Node* term_link;

  char char_to_parent;
  bool is_terminal;
  // Ids of patterns that this terminal vertice (if it is) covers.
  std::vector<unsigned> leaf_pattern_number;
};

class AhoKorasickTrie {
public:
  AhoKorasickTrie(char wildcard_for_patterns, const std::string& temp);
  ~AhoKorasickTrie() = default;

  // Text is read through standart iostream.
  std::vector<unsigned>
  GetNumberOfSubPatternsInPositions();

  unsigned GetNumberOfSubPatterns() const;

private:
  void EndSubPattern(unsigned cur_pattern, Node* cur_node, unsigned cur_ind);

  Node* GetSuffLink(Node* from);
  Node* GetTransition(Node* from, char go_char);
  Node* GetTerminalSuffLink(Node* from);

  std::vector<Node> trie_;
  std::vector<unsigned> pattern_ends_;

  Node* const ROOT;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - -  MAIN - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string patterns;
  std::cin >> patterns;

  AhoKorasickTrie trie('?', patterns);

  const unsigned patterns_count = trie.GetNumberOfSubPatterns();
  std::vector<unsigned> patterns_count_in_pos =
      trie.GetNumberOfSubPatternsInPositions();

  if (patterns_count_in_pos.size() < patterns.length()) {
    return 0;
  }

  for (unsigned i = 0, length = patterns_count_in_pos.size() - patterns.length() + 1;
       i < length; ++i) {
    if (patterns_count_in_pos[i] == patterns_count) {
      std::cout << i << " ";
    }
  }

  return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - DEFINITIONS - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Node::Node()
    : son_node{nullptr}
    , transit{nullptr}
    , parent(nullptr)
    , suff_link(nullptr)
    , term_link(nullptr)
    , char_to_parent(AHO_KORASICK_DEFAULT_CHAR)
    , is_terminal(false)
    , leaf_pattern_number(0) {
  for (unsigned i = 0; i < ALPHABET_SIZE; ++i) {
    son_node[i] = nullptr;
    transit[i] = nullptr;
  }
}

bool Node::IsRoot() const {
  return char_to_parent == AHO_KORASICK_ROOT_CHAR;
}

void Node::MakeRoot() {
  char_to_parent = AHO_KORASICK_ROOT_CHAR;
}

void Node::MakeTerminal(const unsigned pattern_num) {
  is_terminal = true;
  leaf_pattern_number.push_back(pattern_num);
}

void AhoKorasickTrie::EndSubPattern(
    unsigned cur_pattern, Node* cur_node, unsigned cur_ind) {
  cur_node->MakeTerminal(cur_pattern);
  pattern_ends_.push_back(cur_ind);
}

AhoKorasickTrie::AhoKorasickTrie(char wildcard_for_patterns,
                                 const std::string& temp)
    : trie_(0)
    , pattern_ends_(0)
    , ROOT(nullptr) {
  wildcard_for_patterns -= 'a';
  trie_.reserve(temp.size() + 1);
  trie_.emplace_back();
  const_cast<Node*&>(ROOT) = &trie_[0];
  ROOT->MakeRoot();

  unsigned cur_pattern = 0;
  char cur_char = 0;
  bool last_char_was_wildcard = true;
  Node* cur_node = ROOT;
  for (unsigned i = 0, size = temp.length(); i < size; ++i) {
    cur_char = temp[i] - 'a';

    if (cur_char == wildcard_for_patterns) {
      if (!last_char_was_wildcard) {
        EndSubPattern(cur_pattern, cur_node, i - 1);
        ++cur_pattern;
        cur_node = ROOT;
      }

      last_char_was_wildcard = true;
    } else {
      if (cur_node->son_node[cur_char] == nullptr) {
        trie_.emplace_back();

        cur_node->son_node[cur_char] = &trie_.back();

        trie_.back().parent = cur_node;
        trie_.back().char_to_parent = cur_char;
      }

      cur_node = cur_node->son_node[cur_char];

      if (i == size - 1) {
        EndSubPattern(cur_pattern, cur_node, i);
      }

      last_char_was_wildcard = false;
    }
  }
}

std::vector<unsigned>
AhoKorasickTrie::GetNumberOfSubPatternsInPositions() {
  std::vector<unsigned> sub_patterns_count;
  Node* cur_node = ROOT;

  char cur_char = 0;
  for (unsigned i = 0; std::cin >> cur_char; ++i) {
    sub_patterns_count.push_back(0);
    cur_char -= 'a';
    cur_node = GetTransition(cur_node, cur_char);

    for (Node* next_term_node = cur_node;
         !next_term_node->IsRoot();
         next_term_node = GetTerminalSuffLink(next_term_node)) {
      if (next_term_node->is_terminal) {
        for (auto& pattern_id : next_term_node->leaf_pattern_number) {
          if (i >= pattern_ends_[pattern_id]) {
            ++sub_patterns_count[i - pattern_ends_[pattern_id]];
          }
        }
      }
    }
  }

  return sub_patterns_count;
}

unsigned AhoKorasickTrie::GetNumberOfSubPatterns() const {
  return pattern_ends_.size();
}

Node* AhoKorasickTrie::GetSuffLink(Node* from) {
  if (from->suff_link == nullptr) {
    if (from->IsRoot() || from->parent->IsRoot()) {
      from->suff_link = ROOT;
    } else {
      from->suff_link =
          GetTransition(GetSuffLink(from->parent), from->char_to_parent);
    }
  }
  return from->suff_link;
}

Node* AhoKorasickTrie::GetTransition(Node* from, char go_char) {
  if (from->transit[go_char] == nullptr) {
    if (from->son_node[go_char] != nullptr) {
      from->transit[go_char] = from->son_node[go_char];
    } else if (from->IsRoot()) {
      from->transit[go_char] = ROOT;
    } else {
      from->transit[go_char] = GetTransition(GetSuffLink(from), go_char);
    }
  }
  return from->transit[go_char];
}

Node* AhoKorasickTrie::GetTerminalSuffLink(Node* from) {
  if (from->term_link == nullptr) {
    Node* suff_node = GetSuffLink(from); 
    if (suff_node->is_terminal || suff_node->IsRoot()) {
      from->term_link = suff_node;
    } else {
      from->term_link = GetTerminalSuffLink(suff_node);
    }
  }
  return from->term_link;
}
