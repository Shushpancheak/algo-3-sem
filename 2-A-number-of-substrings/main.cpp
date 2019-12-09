/**
 * A. Количество различных подстрок
 * Ограничение времени	0.2 секунды
 * Ограничение памяти	10Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Дана строка длины n. Найти количество ее различных подстрок. Используйте суффиксный массив.
 * Построение суффиксного массива выполняйте за O(n log n). Вычисление количества различных подстрок выполняйте за O(n).
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

// An array that maps every character in a string to some size_t value.
struct StringArray {
  std::vector<size_t> arr;
  std::string string;

  StringArray(const std::string& str);
  virtual ~StringArray() = default;

  StringArray GetInversed() const;

  size_t& operator[](size_t i);
  const size_t& operator[](size_t i) const;
};

class SuffixArray : public StringArray {
public:
  SuffixArray(const std::string& str);
  ~SuffixArray() = default;

private:
  void MakeSuffixArray();
  void SubtractCyclicOffset(const size_t offset);

  std::vector<size_t> StringCountingSort(
      const std::string& str,
      std::vector<size_t>& letters_count);

  void StringCountClassesOfEquivalence(
      const std::string& str,
      const std::vector<size_t>& sorted_positions, 
      std::vector<size_t>& class_eq);

  void ClassesOfEquivalenceCountingSort(
      std::vector<size_t>& class_count,
      std::vector<size_t>& class_eq);

  std::vector<size_t> GetNewClassesOfEquivalence(
      const std::vector<size_t>& current_positions,
      const std::vector<size_t>& old_class_eq,
      const size_t offset);

  const size_t alphabet_size_ = 256;
  size_t total_classes_count = 1;
};

class LcpArray : public StringArray {
public:
  LcpArray(const std::string& str, const SuffixArray& suf_arr);
  ~LcpArray() = default;

private:
  void MakeLcpArray(const SuffixArray& suffix_array);
};


size_t GetNumberOfDifferentSubstrings(const std::string& str);


int main() {
  std::string input_str;
  std::cin >> input_str;

  std::cout << GetNumberOfDifferentSubstrings(input_str);

  return 0;
}


StringArray::StringArray(const std::string& str)
  : arr(str.size(), 0)
  , string(str) {}


size_t& StringArray::operator[](size_t i) {
  return arr[i];
}


const size_t& StringArray::operator[](size_t i) const {
  return arr[i];
}


SuffixArray::SuffixArray(const std::string& str)
  : StringArray(str) {
  MakeSuffixArray();
}


void SuffixArray::SubtractCyclicOffset(const size_t offset) {
  const size_t arr_length = arr.size();

  for (size_t i = 0; i < arr_length; ++i) {
    int element_with_offset = static_cast<int>(arr[i]) -
                              offset;
    if (element_with_offset < 0) {
      element_with_offset += arr_length;
    }
    arr[i] = element_with_offset;
  }
}


std::vector<size_t> SuffixArray::StringCountingSort(
    const std::string& str,
    std::vector<size_t>& letters_count) {
  const size_t str_length = str.length();
  std::vector<size_t> output_sorted_positions(str_length);

  for (size_t i = 0; i < str_length; ++i) {
    ++letters_count[str[i]];
  }
  for (size_t i = 1; i < alphabet_size_; ++i) {
    letters_count[i] += letters_count[i - 1];
  }
  for (size_t i = 0; i < str_length; ++i) {
    output_sorted_positions[--letters_count[str[i]]] = i;
  }

  return output_sorted_positions;
}


void SuffixArray::StringCountClassesOfEquivalence(
    const std::string& str,
    const std::vector<size_t>& sorted_positions, 
    std::vector<size_t>& class_eq) {
  const size_t str_length = str.length();

  class_eq[sorted_positions[0]] = 0;
  for (size_t i = 1; i < str_length; ++i) {
    if (str[sorted_positions[i]] != str[sorted_positions[i - 1]]) {
      ++total_classes_count;
    }
    class_eq[sorted_positions[i]] = total_classes_count - 1;
  }
}


void SuffixArray::ClassesOfEquivalenceCountingSort(
    std::vector<size_t>& class_count,
    std::vector<size_t>& class_eq) {
  const size_t array_size = arr.size();
  std::vector<size_t> output_sorted_positions(array_size);

  std::fill(class_count.begin(), class_count.end(), 0);
  for (size_t i = 0; i < array_size; ++i) {
    ++class_count[class_eq[arr[i]]];
  }
  for (size_t i = 1; i < total_classes_count; ++i) {
    class_count[i] += class_count[i - 1];
  }
  for (int i = static_cast<int>(array_size) - 1; i >= 0; --i) {
    output_sorted_positions[--class_count[class_eq[arr[i]]]] = arr[i];
  }

  arr = output_sorted_positions;
}


std::vector<size_t> SuffixArray::GetNewClassesOfEquivalence(
    const std::vector<size_t>& current_positions,
    const std::vector<size_t>& old_class_eq,
    const size_t offset) {
  const size_t array_size = current_positions.size();
  std::vector<size_t> new_class_eq(old_class_eq.size(), 0);

  new_class_eq[current_positions[0]] = 0;
  total_classes_count = 1;
  for (size_t i = 1; i < array_size; ++i) {
    const size_t cur_second_part = (current_positions[i] + offset) % array_size;
    const size_t prev_second_part = (current_positions[i - 1] + offset) %
                                    array_size;
    if (old_class_eq[current_positions[i]] !=
        old_class_eq[current_positions[i - 1]] ||
        old_class_eq[cur_second_part] != old_class_eq[prev_second_part]) {
      ++total_classes_count;
    }
    new_class_eq[current_positions[i]] = total_classes_count - 1;
  }

  return new_class_eq;
}


void SuffixArray::MakeSuffixArray() {
  const size_t str_length = string.length();
  const size_t max_classes_count = std::max(str_length, alphabet_size_);
  std::vector<size_t> class_count(max_classes_count);
  std::vector<size_t> class_eq(max_classes_count);
  
  // Arr represents current positions of classes of equivalence if they were sorted.
  arr = StringCountingSort(string, class_count);

  // Filling initial classes of equivalence.
  StringCountClassesOfEquivalence(string, arr, class_eq);

  for (size_t offset = 1; offset < str_length; offset *= 2) {
    // For each offset, sort current positions in suffix array according
    // to classes of equivalence and compute new said classes.

    // Using counting sort for pairs of length 2^offset.
    // Since each subarray of length 2^(offset - 1) was previously
    // sorted, we can sort second members of pairs by just subtracting
    // 2^(offset - 1) and the sort first elements.
    SubtractCyclicOffset(offset);
    ClassesOfEquivalenceCountingSort(class_count, class_eq);

    class_eq = GetNewClassesOfEquivalence(arr, class_eq, offset);
  }
}


// Domain should be equal to range.
StringArray StringArray::GetInversed() const {
  const size_t arr_size = arr.size();
  std::vector<size_t> out(arr_size);

  for (size_t i = 0; i < arr_size; ++i) {
    out[arr[i]] = i;
  }

  StringArray out_arr(string);
  out_arr.arr = out;

  return out_arr;
}


LcpArray::LcpArray(const std::string& str, const SuffixArray& suf_arr)
  : StringArray(str) {
  MakeLcpArray(suf_arr);
}


void LcpArray::MakeLcpArray(const SuffixArray& suffix_array) {
  const size_t str_length = string.length();
  // We will need positions of sorted suffixes.
  StringArray pos = suffix_array.GetInversed();

  size_t cur_prefix_len = 0;
  for (size_t i = 0; i < str_length; ++i) {
    if (cur_prefix_len > 0) {
      --cur_prefix_len;
    }

    if (pos[i] == str_length - 1) {
      cur_prefix_len = 0;
    } else {
      const size_t next_suf_start = suffix_array[pos[i] + 1];

      // Comparing current and next in sorted order strings until
      // we reach the end or they are not equal at some character.
      while (std::max(
                 i + cur_prefix_len, next_suf_start + cur_prefix_len) < str_length &&
             string[i + cur_prefix_len] == string[next_suf_start + cur_prefix_len]) {
        ++cur_prefix_len;
      }

      arr[pos[i]] = cur_prefix_len;
    }
  }
}


size_t GetNumberOfDifferentSubstrings(const std::string& str) {
  SuffixArray suf_arr(str);
  LcpArray lcp_arr(str, suf_arr);

  const size_t str_length = str.length();

  size_t substrings_count = 0;

  for (size_t i = 0; i < str_length; ++i) {
    const size_t suffix_length = str_length - suf_arr[i];
    substrings_count += suffix_length - (i == str_length - 1 ? 0 : lcp_arr[i]);
  }

  return substrings_count;
}
