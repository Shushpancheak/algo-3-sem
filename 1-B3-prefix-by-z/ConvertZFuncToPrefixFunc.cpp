#include "ConvertZFuncToPrefixFunc.hpp"

std::vector<size_t> ConvertZFuncToPrefixFunc(const std::vector<size_t>& z_func) {
  std::vector<size_t> prefix_func(z_func.size());

  for (size_t i = 0, z_length = z_func.size(); i < z_length; ++i) {
    for (int j = z_func[i] - 1; j >= 0; --j) {
      if (prefix_func[i + j] > 0) {
        break;
      } else {
        prefix_func[i + j] = j + 1;
      }
    }
  }

  return prefix_func;
}
