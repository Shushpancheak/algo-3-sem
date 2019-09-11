/**
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

std::vector<size_t> GetAllSubstringPositions(const std::string& str, const std::string& substr);

int main() {
  std::string str;
  std::string substr;

  std::cin >> substr;
  std::cin >> str;

  for (auto& ind : GetAllSubstringPositions(str, substr)) {
    std::cout << ind << " ";
  }

  return 0;
}

std::vector<size_t> GetAllSubstringPositions(const std::string& str, const std::string& substr) {
  std::vector<size_t> res;
  std::vector<size_t> prefix_func(substr.length());

  // Calculating prefix function for substring.
  for (size_t i = 1; i < substr.length(); ++i) {
    size_t cur_prefix_size = prefix_func[i - 1];
    // "Jumping"
    while (cur_prefix_size) {
      if (substr[i] == substr[cur_prefix_size]) {
        prefix_func[i] = cur_prefix_size + 1;
        break;
      }
      cur_prefix_size = prefix_func[cur_prefix_size - 1];
    }
    if (!prefix_func[i]) {
      prefix_func[i] = substr[0] == substr[i] ? 1 : 0;
    }
  }

  // Calculating prefix size for every character in string.
  size_t cur_prefix_size = 0;
  for (size_t i = 0; i < str.length(); ++i) {
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
      res.push_back(i - substr.length() + 1);
    }
  }

  return res;
}
