/*
 * B3. Перевести Z-функцию в префикс-функцию.
 *
 * T(n) = O(n);
 * M(n) = O(n);
 */
#include <vector>
#include <iostream>
#include "ConvertZFuncToPrefixFunc.hpp"

int main() {
  std::vector<size_t> z_func;

  std::cout << "Please input the Z-function sequence: ";
  size_t next_elem;
  while (std::cin >> next_elem) {
    z_func.push_back(next_elem);
  }
  std::cout << "\n";

  for (auto& elem : ConvertZFuncToPrefixFunc(z_func)) {
    std::cout << elem << " ";
  }
  return 0;
}
