/**
 * B1. Строка по префикс-функции
 * Ограничение времени	0.25 секунд
 * Ограничение памяти	64Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Найти лексикографически-минимальную строку, построенную по префикс-функции, в алфавите a-z.
 *
 *  T(n) = O(n);
 *  M(n) = O(n).
 */
#include <iostream>
#include <vector>
#include <string>
#include <set>

char GetMinAvailableSymbol(const std::set<char>& forbidden_symbols);
std::string GetMinStrByPrefixFunc(const std::vector<unsigned>& prefix_func);

int main() {
  std::vector<unsigned> prefix_func;

  unsigned elem = 0;
  while (std::cin >> elem) {
    prefix_func.push_back(elem);
  }

  std::cout << GetMinStrByPrefixFunc(prefix_func);

  return 0;
}

char GetMinAvailableSymbol(const std::set<char>& forbidden_symbols) {
  char next_symbol = 'b';
  for (; ; ++next_symbol) {
    if (forbidden_symbols.find(next_symbol) == forbidden_symbols.end()) {
      break;
    }
  }
  return next_symbol;
}

std::string GetMinStrByPrefixFunc(const std::vector<unsigned>& prefix_func) {
  std::string min_str;
  min_str.reserve(prefix_func.size());

  if (!prefix_func.empty()) {
    min_str = "a";
  }
  for (size_t i = 1; i < prefix_func.size(); ++i) {
    if (prefix_func[i]) {
      min_str += min_str[prefix_func[i] - 1];
    } else {
      std::set<char> forbidden_symbols;

      size_t cur_prefix_size = prefix_func[i - 1];
      while (cur_prefix_size) {
        forbidden_symbols.insert(min_str[cur_prefix_size]);
        cur_prefix_size = prefix_func[cur_prefix_size - 1];
      }

      min_str += GetMinAvailableSymbol(forbidden_symbols);
    }
  }

  return min_str;
}
