﻿/**
 * B2. Строка по Z-функции
 * Ограничение времени	0.15 секунд
 * Ограничение памяти	  64Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Найти лексикографически-минимальную строку, построенную по z-функции, в алфавите a-z.
 *
 *  T(n) = O(n);
 *  M(n) = O(n).
 */
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

std::string GetMinStrByZFunc(const std::vector<int>& z_func);
char GetMinAvailableSymbol(const std::string& str, const std::vector<size_t>& forbidden_pos);

int main() {
  std::vector<int> prefix_func;

  int elem = 0;
  while (std::cin >> elem) {
    prefix_func.push_back(elem);
  }

  std::cout << GetMinStrByZFunc(prefix_func);

  return 0;
}

char GetMinAvailableSymbol(const std::string& str, const std::vector<size_t>& forbidden_pos) {
  std::set<char> forbidden_symbols;
  for (auto& pos : forbidden_pos) {
    forbidden_symbols.insert(str[pos]);
  }
  char next_symbol = 'b';
  for (; ; ++next_symbol) {
    if (forbidden_symbols.find(next_symbol) == forbidden_symbols.end()) {
      break;
    }
  }
  return next_symbol;
}

std::string GetMinStrByZFunc(const std::vector<int>& z_func) {
  std::string min_str;
  min_str.reserve(z_func.size());

  if (!z_func.empty()) {
    min_str = "a";
  }

  // Borders of the right-most window
  int right = -1, left = 0;

  // Positions of symbols which must not be written once we
  // stumble upon 0's. They are the next symbols of prefixes
  // in current window. If we would've written one of these symbols,
  // previous z-function values would've been incorrect.
  std::vector<size_t> forbidden_pos;

  // Using i as int because we could have negative borders by default.
  for (int i = 1; i < static_cast<int>(z_func.size()); ++i) {
    if (z_func[i]) {
      // We are standing at the start of prefix, so adding 'a'.
      min_str += 'a';

      // If the current window is smaller than this, create the new
      // window and clear forbidden pos.
      if (z_func[i] + i - 1 > right) {
        left = i;
        right = i + z_func[i] - 1;
        forbidden_pos.clear();
      }

      // Add to forbidden_pos only if the current window touches
      // the right border of our window.
      if (z_func[i] + i - 1 == right) {
        forbidden_pos.push_back(z_func[i]);
      }
    } else {
      // If we are within the current window, just copy...
      if (i <= right) {
        min_str += min_str[i - left];
      } else {
        // Choosing the next symbol that won't break
        // previous z-func values.
        min_str += GetMinAvailableSymbol(min_str, forbidden_pos);

        // Creating new window with 0 size.
        left = right = i;
        forbidden_pos.clear();
      }
    }
  }

  return min_str;
}
