﻿/**
 * A. Поиск подстроки
 * Ограничение времени	0.1 секунда
 * Ограничение памяти	32Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Найдите все вхождения шаблона в строку. Длина шаблона – p, длина строки – n. Время O(n + p), доп. память – O(p).
 * p <= 30000, n <= 300000.
 * Использовать один из методов:
 * - С помощью префикс-функции;
 * - С помощью z-функции.
 *
 *  T(n, p) = O(n + p)
 *  M(n, p) = O(p).
 */
#include <iostream>
#include <vector>
#include <string>

std::vector<size_t> CalculatePrefixFunc(const std::string& str);
std::vector<size_t> GetAllSubstringPositions(const std::string& str, const std::string& substr);

int main() {
  std::string str;
  std::string substr;

  std::cin >> substr;
  std::cin >> str;

  for (auto& pos : GetAllSubstringPositions(str, substr)) {
    std::cout << pos << " ";
  }

  return 0;
}

std::vector<size_t> CalculatePrefixFunc(const std::string& str) {
  std::vector<size_t> prefix_func(str.length());
  for (size_t i = 1, str_length = str.length(); i < str_length(); ++i) {
    size_t cur_prefix_size = prefix_func[i - 1];
    // "Jumping"
    while (cur_prefix_size) {
      if (str[i] == str[cur_prefix_size]) {
        prefix_func[i] = cur_prefix_size + 1;
        break;
      }
      cur_prefix_size = prefix_func[cur_prefix_size - 1];
    }
    if (!prefix_func[i]) {
      prefix_func[i] = str[0] == str[i] ? 1 : 0;
    }
  }
  return prefix_func;
}

std::vector<size_t> GetAllSubstringPositions(const std::string& str, const std::string& substr) {
  std::vector<size_t> substr_pos;

  // Calculating prefix function for substring.
  std::vector<size_t> prefix_func = CalculatePrefixFunc(substr);

  // Calculating prefix size for every character in string.
  size_t cur_prefix_size = 0;
  for (size_t i = 0, str_length = str.length(); i < str_length; ++i) {
    // "Jumping"
    while (cur_prefix_size) {
      if (cur_prefix_size < substr.length() && str[i] == substr[cur_prefix_size]) {
        ++cur_prefix_size;
        break;
      }
      cur_prefix_size = prefix_func[cur_prefix_size - 1];
    }
    if (!cur_prefix_size) {
      cur_prefix_size = substr[0] == str[i] ? 1 : 0;
    }
    if (cur_prefix_size == substr.length()) {
      substr_pos.push_back(i - substr.length() + 1);
    }
  }

  return substr_pos;
}
