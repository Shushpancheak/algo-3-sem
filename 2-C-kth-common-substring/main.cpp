#include <iostream>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

std::vector<size_t> GetSuffixArray(const std::string& str);
std::vector<size_t> GetLcpArray(const std::string& str,
                                const std::vector<size_t>& suffix_array);
std::string GetKthCommonSubstring(const std::string& str1,
                                  const std::string& str2,
                                  size_t k);


int main() {
  std::string input_str_1, input_str_2;
  size_t k;
  std::cin >> input_str_1 >> input_str_2 >> k;

  std::cout << GetKthCommonSubstring(input_str_1, input_str_2, k);

  return 0;
}


std::vector<size_t> SubtractCyclicOffset(std::vector<size_t>& arr,
                                         const size_t offset) {
  const size_t arr_length = arr.size();
  std::vector<size_t> new_arr(arr_length);
  for (size_t i = 0; i < arr_length; ++i) {
    int element_with_offset = static_cast<int>(arr[i]) -
                              offset;
    if (element_with_offset < 0) {
      element_with_offset += arr_length;
    }
    new_arr[i] = element_with_offset;
  }
  return new_arr;
}


std::vector<size_t> StringCountingSort(
    const std::string& str,
    const size_t alphabet_size,
    std::vector<size_t>& letters_count) {
  const size_t str_length = str.length();
  std::vector<size_t> output_sorted_positions(str_length);

  for (size_t i = 0; i < str_length; ++i) {
    ++letters_count[str[i]];
  }
  for (size_t i = 1; i < alphabet_size; ++i) {
    letters_count[i] += letters_count[i - 1];
  }
  for (size_t i = 0; i < str_length; ++i) {
    output_sorted_positions[--letters_count[str[i]]] = i;
  }

  return output_sorted_positions;
}


void StringCountClassesOfEquivalence(
    const std::string& str,
    const std::vector<size_t>& sorted_positions, 
    std::vector<size_t>& class_eq,
    size_t& total_classes_count) {
  const size_t str_length = str.length();

  class_eq[sorted_positions[0]] = 0;
  for (size_t i = 1; i < str_length; ++i) {
    if (str[sorted_positions[i]] != str[sorted_positions[i - 1]]) {
      ++total_classes_count;
    }
    class_eq[sorted_positions[i]] = total_classes_count - 1;
  }
}


std::vector<size_t> ClassesOfEquivalenceCountingSort(
    std::vector<size_t>& class_count,
    std::vector<size_t>& class_eq,
    std::vector<size_t>& old_positions,
    size_t total_classes_count) {
  const size_t array_size = old_positions.size();
  std::vector<size_t> output_sorted_positions(array_size);

  std::fill(class_count.begin(), class_count.end(), 0);
  for (size_t i = 0; i < array_size; ++i) {
    ++class_count[class_eq[old_positions[i]]];
  }
  for (size_t i = 1; i < total_classes_count; ++i) {
    class_count[i] += class_count[i - 1];
  }
  for (int i = static_cast<int>(array_size) - 1; i >= 0; --i) {
    output_sorted_positions[--class_count[class_eq[old_positions[i]]]] = old_positions[i];
  }

  return output_sorted_positions;
}


std::vector<size_t> GetNewClassesOfEquivalence(
    const std::vector<size_t>& current_positions,
    const std::vector<size_t>& old_class_eq,
    size_t& total_classes_count,
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


std::vector<size_t> GetSuffixArray(const std::string& str) {
  const size_t alphabet_size = 256;
  const size_t str_length = str.length();
  const size_t max_classes_count = std::max(str_length, alphabet_size);
  // Current positions of classes of equivalence if they were sorted.
  std::vector<size_t> current_positions(str_length);
  std::vector<size_t> class_count(max_classes_count);
  std::vector<size_t> class_eq(max_classes_count);
  std::vector<size_t> new_positions(str_length);
  size_t total_classes_count = 1;

  current_positions = StringCountingSort(str, alphabet_size, class_count);

  StringCountClassesOfEquivalence(str, current_positions,
                                  class_eq, total_classes_count);

  
  for (size_t offset = 1; offset < str_length; offset *= 2) {
    new_positions = SubtractCyclicOffset(current_positions, offset);
    current_positions = ClassesOfEquivalenceCountingSort(class_count, class_eq,
                                     new_positions, total_classes_count);
    
    class_eq = GetNewClassesOfEquivalence(current_positions, class_eq, 
                                          total_classes_count, offset);
  }

  return current_positions;
}


// Domain should be equal to range.
std::vector<size_t> GetInversed(const std::vector<size_t>& arr) {
  const size_t arr_size = arr.size();
  std::vector<size_t> out(arr_size);

  for (size_t i = 0; i < arr_size; ++i) {
    out[arr[i]] = i;
  }

  return out;
}


std::vector<size_t> GetLcpArray(const std::string& str,
                                const std::vector<size_t>& suffix_array) {
  const size_t str_length = str.length();
  std::vector<size_t> lcp(str_length);
  std::vector<size_t> pos = GetInversed(suffix_array);

  size_t cur_prefix = 0;
  for (size_t i = 0; i < str_length; ++i) {
    if (cur_prefix > 0) {
      --cur_prefix;
    }
    if (pos[i] == str_length - 1) {
      cur_prefix = 0;
    } else {
      const size_t next_suf = suffix_array[pos[i] + 1];
      while (std::max(i + cur_prefix, next_suf + cur_prefix) < str_length &&
             str[i + cur_prefix] == str[next_suf + cur_prefix]) {
        ++cur_prefix;
      }
      lcp[pos[i]] = cur_prefix;
    }
  }

  return lcp;
}


std::string GetKthCommonSubstring(const std::string& str1,
                                  const std::string& str2,
                                  size_t k) {
  const std::string str = str1 + '\1' + str2 + '\2';

  std::vector<size_t> suf_arr = GetSuffixArray(str);
  std::vector<size_t> lcp_arr = GetLcpArray(str, suf_arr);

  const size_t str_length = str.length();
  const size_t str2_start = str1.length() + 1;
  size_t substrings_count = 0;
  size_t last_used_lcp = 0;

  for (size_t i = 0; i < str_length - 1; ++i) {
    if ((suf_arr[i] < str2_start && suf_arr[i + 1] >= str2_start) ||
        (suf_arr[i] >= str2_start && suf_arr[i + 1] < str2_start)) {
      const size_t plus_diff = lcp_arr[i] - std::min(lcp_arr[i], last_used_lcp);
      substrings_count += plus_diff;
      if (substrings_count >= k) {
        const size_t old_substrings_count = substrings_count - plus_diff;
        const size_t out_size = k - old_substrings_count + last_used_lcp;

        return str.substr(suf_arr[i], out_size);
      }

      last_used_lcp = lcp_arr[i];
    } else if (lcp_arr[i] < last_used_lcp) {
      last_used_lcp = lcp_arr[i];
    }
  }

  return "-1";
}
