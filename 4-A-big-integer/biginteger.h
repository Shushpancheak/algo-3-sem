#ifndef BIG_INTEGER
#define BIG_INTEGER

#include <iostream>
#include <vector>
#include <string>

// To use ssize_t in MVC.
#ifdef _MSC_VER
#include <basetsd.h>
typedef SSIZE_T ssize_t;
#endif

static size_t counter = 0;
static const size_t CNT = 160;

namespace details {

// The type which size is enough for storing +-2 * digit
// in the maximum base 10^n for this type (which will be calculated later).
#define DIGIT_TYPE     int32_t
// A type that is capable of holding digit * digit.
#define MUL_DIGIT_TYPE int64_t
// Since c++11 compiler doesn't support proper constexprs, i'm just hard coding
// the maximum base for the int32_t.
#define DIGIT_BASE     1e9

struct QuotientAndRemainder;

// Fill trailing zeros so that the total length will be len
void FillTrailingZeros(std::vector<DIGIT_TYPE>& vec, size_t len);
void FillTrailingZeros(std::string& str,             size_t len);
void FillFirstZeros(std::vector<DIGIT_TYPE>& vec, size_t len);
void FillFirstZeros(std::string& str,             size_t len);
void DeleteTrailingZeros(std::vector<DIGIT_TYPE>& vec);

}


class BigInteger {
public:
  BigInteger();
  BigInteger(int integer);
  BigInteger(const BigInteger& bigint) = default;
  BigInteger(BigInteger&& bigint) noexcept;
  ~BigInteger() = default;

  BigInteger  operator- ();
  BigInteger& operator= (const BigInteger&  bigint) = default;
  BigInteger& operator= (      BigInteger&& bigint) noexcept;
  BigInteger& operator+=(const BigInteger&  bigint);
  BigInteger& operator-=(const BigInteger&  bigint);
  BigInteger& operator/=(const BigInteger&  bigint);
  BigInteger& operator*=(const BigInteger&  bigint);
  BigInteger& operator%=(const BigInteger&  bigint);
  BigInteger& operator++();
  BigInteger  operator++(int);
  BigInteger& operator--();
  BigInteger  operator--(int);

  // ReSharper disable once CppInconsistentNaming
  std::string toString() const;

  explicit operator int() const;
  explicit operator bool() const;

  friend BigInteger operator+(const BigInteger& val_1, const BigInteger& val_2);
  friend BigInteger operator-(const BigInteger& val_1, const BigInteger& val_2);
  friend BigInteger operator/(const BigInteger& val_1, const BigInteger& val_2);
  friend BigInteger operator*(const BigInteger& val_1, const BigInteger& val_2);
  friend BigInteger operator%(const BigInteger& val_1, const BigInteger& val_2);

  friend bool operator==(const BigInteger& val_1, const BigInteger& val_2);
  friend bool operator!=(const BigInteger& val_1, const BigInteger& val_2);
  friend bool operator> (const BigInteger& val_1, const BigInteger& val_2);
  friend bool operator< (const BigInteger& val_1, const BigInteger& val_2);
  friend bool operator>=(const BigInteger& val_1, const BigInteger& val_2);
  friend bool operator<=(const BigInteger& val_1, const BigInteger& val_2);

  friend std::istream&  operator>>(std::istream&  stream,       BigInteger& bigint);
  friend std::istream&& operator>>(std::istream&& stream,       BigInteger& bigint);
  friend std::ostream&  operator<<(std::ostream&  stream, const BigInteger& bigint);
  friend std::ostream&& operator<<(std::ostream&& stream, const BigInteger& bigint);

private:
  static const DIGIT_TYPE BASE = DIGIT_BASE;

  // Used both in operator+= and -=, since we don't want to
  // copy-paste code and to copy object (with unary -) to pass
  // to operator+= in operator-=..
  void Add(const BigInteger& bigint, bool is_bigint_negative);

  // Used in division and module to reduce copy-pasting...
  friend details::QuotientAndRemainder Divide(const BigInteger& dividend,
                                              const BigInteger& divisor);
  // Finds a digit such that dividend - digit * divisor >= 0, but dividend -
  // - (digit + 1) * divisor < 0 via bin search.
  friend DIGIT_TYPE FindMinimalDigitForDivision(const BigInteger& dividend,
                                                const BigInteger& divisor);

  // Split number in two
  std::pair<BigInteger, BigInteger> SplitNumber(const BigInteger& bigint,
                                                size_t digit_ind);
  // Does not affect the sign.
  BigInteger KaratsubaMultiply(const BigInteger& bigint_1, const BigInteger& bigint_2);

  bool is_negative_;
  // Starting from least significant digit.
  std::vector<DIGIT_TYPE> digits_;
};

struct details::QuotientAndRemainder {
  BigInteger quotient;
  BigInteger remainder;
};


template<typename T>
static T Min(const T& first, const T& second) {
  return first < second ? first : second;
}


template<typename T>
static T Max(const T& first, const T& second) {
  return first > second ? first : second;
}


inline BigInteger operator+(const BigInteger& val_1, const BigInteger& val_2) {
  BigInteger bigint(val_1);
  return bigint += val_2;
}


inline BigInteger operator-(const BigInteger& val_1, const BigInteger& val_2) {
  BigInteger bigint(val_1);
  return bigint -= val_2;
}


inline BigInteger operator/(const BigInteger& val_1, const BigInteger& val_2) {
  BigInteger bigint(val_1);
  return bigint /= val_2;
}


inline BigInteger operator*(const BigInteger& val_1, const BigInteger& val_2) {
  BigInteger bigint(val_1);
  return bigint *= val_2;
}


inline BigInteger operator%(const BigInteger& val_1, const BigInteger& val_2) {
  BigInteger bigint(val_1);
  return bigint %= val_2;
}


inline bool operator==(const BigInteger& val_1, const BigInteger& val_2) {
  // Both zeros.
  if (!val_1 && !val_2) {
    return true;
  }

  return val_1.is_negative_ == val_2.is_negative_ && val_1.digits_ == val_2.digits_;
}


inline bool operator!=(const BigInteger& val_1, const BigInteger& val_2) {
  return !(val_1 == val_2);
}


inline bool operator>(const BigInteger& val_1, const BigInteger& val_2) {
  if (val_1.digits_.size() > val_2.digits_.size()) {
    if (!val_1.is_negative_) {
      return true;
    } else {
      return false;
    }
  } else if (val_1.digits_.size() < val_2.digits_.size()) {
    if (!val_2.is_negative_) {
      return false;
    } else {
      return true;
    }
  }

  // If both are zero.
  if (!val_1) {
    return false;
  }

  // length_1 == length_2 now.
  if (!val_1.is_negative_ && val_2.is_negative_) {
    return true;
  } else if (val_2.is_negative_ && !val_1.is_negative_) {
    return false;
  }

  // length_1 == length_2 and sign_1 == sign_2 now.
  for (ssize_t i = val_1.digits_.size() - 1; i >= 0; --i) {
    if (val_1.digits_[i] == val_2.digits_[i]) {
      continue;
    }

    const bool first_is_bigger_than_second = val_1.digits_[i] > val_2.digits_[i];
    // v_1 > v_2 <=> -v_1 < -v_2.
    if (val_1.is_negative_ ? !first_is_bigger_than_second :
                              first_is_bigger_than_second){
      return true;
    } else {
      return false;
    }
  }

  return false;
}


inline bool operator<(const BigInteger& val_1, const BigInteger& val_2) {
  return !(val_1 >= val_2);
}


inline bool operator>=(const BigInteger& val_1, const BigInteger& val_2) {
  return val_1 == val_2 || val_1 > val_2;
}


inline bool operator<=(const BigInteger& val_1, const BigInteger& val_2) {
  return !(val_1 > val_2);
}


inline std::istream& operator>>(std::istream& stream, BigInteger& bigint) {
  bigint = 0;

  std::string num_str;
  stream >> num_str;

  bigint.is_negative_ = num_str[0] == '-';

  DIGIT_TYPE    cur_digit  = 0;
  DIGIT_TYPE    cur_pow_10 = 1;
  const ssize_t end        = (bigint.is_negative_ ? 1 : 0);
  for (ssize_t i = num_str.length() - 1; i >= end; --i) {
    if (cur_pow_10 == BigInteger::BASE) {
      bigint.digits_.push_back(cur_digit);

      cur_digit  = num_str[i] - '0';
      cur_pow_10 = 10;
    } else {
      cur_digit  += (num_str[i] - '0') * cur_pow_10;
      cur_pow_10 *= 10;
    }
  }

  bigint.digits_.push_back(cur_digit);

  return stream;
}


inline std::istream&& operator>>(std::istream&& stream, BigInteger& bigint) {
  return std::move(stream >> bigint);
}


inline std::ostream& operator<<(std::ostream& stream, const BigInteger& bigint) {
  return stream << bigint.toString();
}


inline std::ostream&& operator<<(std::ostream&& stream, const BigInteger& bigint) {
  return std::move(stream << bigint);
}


inline void details::FillTrailingZeros(std::vector<DIGIT_TYPE>& vec, size_t len) {
  if (len <= vec.size()) {
    return;
  }

  vec.resize(len, 0);
}


inline void details::FillTrailingZeros(std::string& str, size_t len) {
  if (len <= str.size()) {
    return;
  }

  str.resize(len, '0');
}


inline void details::FillFirstZeros(std::vector<DIGIT_TYPE>& vec, size_t len) {
  vec.insert(vec.begin(), len - vec.size(), 0);
}


inline void details::FillFirstZeros(std::string& str, size_t len) {
  str.insert(str.begin(), len - str.size(), '0');
}


inline void details::DeleteTrailingZeros(std::vector<DIGIT_TYPE>& vec) {
  for (ssize_t i = vec.size() - 1; i >= 0; --i) {
    if (vec[i] == 0) {
      vec.pop_back();
    } else {
      break;
    }
  }
}


inline DIGIT_TYPE FindMinimalDigitForDivision(
    const BigInteger& dividend,
    const BigInteger& divisor) {
  // O(n log BASE).

  DIGIT_TYPE lower = 0, upper = BigInteger::BASE;
  DIGIT_TYPE middle = 0;
  BigInteger cur_bigint = 0;

  while (upper - lower > 1) {
    middle = (upper + lower) / 2;
    cur_bigint = divisor * static_cast<BigInteger>(middle);

    if (cur_bigint > dividend) {
      upper = middle;
    } else {
      lower = middle;
    }
  }

  return lower;
}


inline BigInteger::BigInteger()
  : is_negative_(false)
  , digits_(0) {}


inline BigInteger::BigInteger(int integer)
  : is_negative_(integer < 0)
  , digits_(0) {
  int64_t integer64 = integer; // Needed for case integer == -MAX_INT - 1,
                               // in which case abs(integer) ?

  if (integer64 < 0) {
    integer64 *= -1;
  }

  while (integer64 != 0) {
    digits_.push_back(integer64 % BASE);
    integer64 /= BASE;
  }
}


inline BigInteger::BigInteger(BigInteger&& bigint) noexcept
  : is_negative_(bigint.is_negative_)
  , digits_(std::move(bigint.digits_)) {
  bigint.digits_      = std::vector<DIGIT_TYPE>();
  bigint.is_negative_ = false;
}


inline BigInteger BigInteger::operator-() {
  BigInteger new_int(*this);
  new_int.is_negative_ = !is_negative_;

  return new_int;
}


inline BigInteger& BigInteger::operator=(BigInteger&& bigint) noexcept {
  if (this == &bigint) {
    return *this;
  }

  digits_      = std::move(bigint.digits_);
  is_negative_ = bigint.is_negative_;

  bigint.digits_      = std::vector<DIGIT_TYPE>();
  bigint.is_negative_ = false;

  return *this;
}


inline BigInteger& BigInteger::operator+=(const BigInteger& bigint) {
  Add(bigint, bigint.is_negative_);
  return *this;
}


inline BigInteger& BigInteger::operator-=(const BigInteger& bigint) {
  Add(bigint, !bigint.is_negative_);
  return *this;
}


inline BigInteger& BigInteger::operator/=(const BigInteger& bigint) {
  return *this = Divide(*this, bigint).quotient;
}


inline BigInteger& BigInteger::operator*=(const BigInteger& bigint) {
  const bool this_is_negative = this->is_negative_;

  // One of operands is zero.
  if (!(*this) || !bigint) {
    digits_ = std::vector<DIGIT_TYPE>(); // Zeroing this value.
    is_negative_ = false;
    return *this;
  }

  const size_t len_1 = digits_.size();
  const size_t len_2 = bigint.digits_.size();

  // If one of operands has size of one, then do a trivial O(n) multiplication
  // by a single digit.
  if (len_1 == 1 || len_2 == 1) {
    const bool       this_is_digit = len_1 == 1;
    const DIGIT_TYPE digit         = (this_is_digit ? digits_[0] : bigint.digits_[0]);
    const size_t     len           = (this_is_digit ? len_2      : len_1);
    MUL_DIGIT_TYPE   carry         = 0;
    MUL_DIGIT_TYPE   cur_digit     = 0;

    details::FillTrailingZeros(digits_, len);

    for (size_t i = 0; i < len; ++i) {
      cur_digit  = static_cast<MUL_DIGIT_TYPE>(digit) *
                   (this_is_digit ? bigint.digits_[i] : digits_[i]) + carry;
      carry      = cur_digit / BASE;
      digits_[i] = cur_digit % BASE;
    }

    if (carry > 0) {
      digits_.push_back(carry);
    }
  } else {
    *this = KaratsubaMultiply(*this, bigint);
  }

  // Compute the resulting sign.
  is_negative_ = this_is_negative ^ bigint.is_negative_;

  return *this;
}


inline BigInteger& BigInteger::operator%=(const BigInteger& bigint) {
  return *this = Divide(*this, bigint).remainder;
}


inline BigInteger& BigInteger::operator++() {
  return *this += 1;
}


inline BigInteger BigInteger::operator++(int) {
  BigInteger old_bigint(*this);
  ++*this;
  return old_bigint;
}


inline BigInteger& BigInteger::operator--() {
  return *this -= 1;
}


inline BigInteger BigInteger::operator--(int) {
  BigInteger old_bigint(*this);
  --*this;
  return old_bigint;
}


inline std::string BigInteger::toString() const {
  const size_t length                  = digits_.size();
  const size_t digit_length_in_base_10 = log10(BASE);

  if (length == 0) {
    return "0";
  }

  std::string res_str;
  res_str.reserve(length * digit_length_in_base_10 + 1);
  if (is_negative_) {
    res_str += "-";
  }

  std::string cur_str = std::to_string(digits_[length - 1]);
  res_str += cur_str;
  for (ssize_t i = length - 2; i >= 0; --i) {
    cur_str = std::to_string(digits_[i]);
    details::FillFirstZeros(cur_str, digit_length_in_base_10);
    res_str += cur_str;
  }

  return res_str;
}


inline BigInteger::operator int() const {
  int cur_pow = 1;
  int res     = 0;

  for (auto& digit : digits_) {
    res     += cur_pow * digit;
    cur_pow *= BASE;
  }

  return (is_negative_ ? -1 : 1) * res;
}


inline BigInteger::operator bool() const {
  for (size_t i = 0, len = digits_.size(); i < len; ++i) {
    if (digits_[i] != 0) {
      return true;
    }
  }

  return false;
}


inline void BigInteger::Add(const BigInteger& bigint, const bool is_bigint_negative) {
  const bool   same_sign = is_negative_ == is_bigint_negative;
  DIGIT_TYPE   carry     = 0;
  DIGIT_TYPE   cur_digit = 0;
  const size_t len_2     = bigint.digits_.size();

  // Filling with zeros so that len_1 >= len_2.
  details::FillTrailingZeros(digits_, len_2);
  const size_t len_1 = digits_.size();

  // Just a column addition.
  for (size_t i = 0; i < len_1; ++i) {
    if (i >= len_2) {
      cur_digit = digits_[i] + carry;
    } else {
      cur_digit = digits_[i] + bigint.digits_[i] * (same_sign ? 1 : -1) + carry;
    }

    if (cur_digit >= BASE) {
      carry = 1;
      cur_digit %= BASE;
    } else if (cur_digit < 0) {
      carry = -1;
      cur_digit = BASE + cur_digit;
    } else {
      carry = 0;
    }

    digits_[i] = cur_digit;

    if (i >= len_2 && carry == 0) {
      // No further changes will be done, so just break.
      break;
    }
  }

  // At this point, we either reached len_1 or left earlier but with carry == 0.
  if (carry == 1) {
    digits_.push_back(1);
  } else if (carry == -1) {
    // We wanted to compute a - b, but if carry is -1 in the end,
    // that means we "added" 10^n to a and wanted to carry the 1 away,
    // but that will not be a - b. Instead, consider d = 10^n + a - b.
    // Then, a - b = -(10^n - d).
    // So we will make such thing with given d:
    //  1  0  0  ... 0           
    // -                         
    //     d1 d2 ... d_{n - 1}   
    // --------------------------
    // 0 9-d1 9-d2...10-d_{n - 1}

    is_negative_ = !is_negative_;
    for (size_t i = 0; i < len_1; ++i) {
      digits_[i] = BASE - digits_[i] - (i == 0 ? 0 : 1);
    }
  }

  details::DeleteTrailingZeros(digits_);
}


inline std::pair<BigInteger, BigInteger> BigInteger::SplitNumber(
    const BigInteger& bigint,
    size_t digit_ind) {
  // O(n).

  BigInteger first  = 0;
  BigInteger second = 0;

  const std::vector<DIGIT_TYPE> first_vec (bigint.digits_.begin(),
                                           bigint.digits_.begin() + digit_ind);
  const std::vector<DIGIT_TYPE> second_vec(bigint.digits_.begin() + digit_ind,
                                           bigint.digits_.end());

  first. digits_ = first_vec;
  second.digits_ = second_vec;

  return std::make_pair(first, second);
}


inline BigInteger BigInteger::KaratsubaMultiply(
    const BigInteger& bigint_1,
    const BigInteger& bigint_2) {
  const size_t len_1 = bigint_1.digits_.size();
  const size_t len_2 = bigint_2.digits_.size();

  if (len_1 <= 1 || len_2 <= 1) {
    // One-digit multiplication is defined there (O(n)).
    return bigint_1 * bigint_2;
  }

  const size_t middle = Min(len_1, len_2) / 2;

  const std::pair<BigInteger, BigInteger> parts_1 = SplitNumber(bigint_1, middle);
  const std::pair<BigInteger, BigInteger> parts_2 = SplitNumber(bigint_2, middle);

  // bigint_1 = upp_1 * BASE^m + low_1,
  // bigint_2 = upp_2 * BASE^m + low_2,
  // bigint_1 * bigint_2 = (upp_1 * BASE^m + low_1) * (upp_2 * BASE^m + low_2) =
  //                     = first * B^(2m) + second * B^m + third.
  // first  = upp_1 * upp_2,
  // third  = low_1 * low_2,
  // second = low_1 * upp_2 + low_2 * upp_1 =
  //        = (low_1 + upp_1)(low_2 + upp_2) - first - second.
  BigInteger first  = KaratsubaMultiply(parts_1.second, parts_2.second);
  BigInteger third  = KaratsubaMultiply(parts_1.first, parts_2.first);
  BigInteger second = KaratsubaMultiply(parts_1.first + parts_1.second,
                                        parts_2.first + parts_2.second) - first
                                                                        - third;

  details::FillFirstZeros(first.digits_,  first.digits_. size() + middle * 2);
  details::FillFirstZeros(second.digits_, second.digits_.size() + middle);

  return first + second + third;
}


inline details::QuotientAndRemainder Divide(
    const BigInteger& dividend,
    const BigInteger& divisor) {
  // Dividing using the school method with asymptotic O(n^2):
  //
  // 1. Shift the divisor so that it was aligned with dividend        = O(n);
  // 2. Find the minimal number such that dividend - number * divisor
  //    >= 0, but dividend - (number + 1) * divisor < 0 and push
  //    this number (which should be a digit in the given base) to
  //    quotient bigint. Since the multiplication by a digit has only
  //    O(n) complexity and we will find that number with bin search,
  //    the resulting asymptotic                                      = O(n log BASE)
  // 3. Subtract that number from dividend and go to 1.               = O(n)
  //
  // Doing such steps no more than n times, we get the resulting
  // asymptotic of O(n^2 log BASE).

  details::QuotientAndRemainder result = {0, dividend};
  BigInteger cur_divisor       = divisor;
  BigInteger divisor_abs_value = divisor;
  divisor_abs_value.is_negative_ = false;

  if (!dividend) {
    return result;
  } else if (!cur_divisor) {
    throw(std::invalid_argument("Division by zero in BigInteger."));
  }

  // Compute the resulting sign.
  const bool is_quotient_negative  = dividend.is_negative_ ^ cur_divisor.is_negative_;
  const bool is_remainder_negative = dividend.is_negative_;

  const size_t dividend_length = dividend.   digits_.size();
  const size_t divisor_length  = cur_divisor.digits_.size();

  if (dividend_length < divisor_length) {
    return result;
  }

  // Setting size of the quotient.
  const size_t dividend_divisor_len_diff = dividend_length - divisor_length;
  details::FillTrailingZeros(result.quotient.digits_, dividend_divisor_len_diff + 1);
  ssize_t cur_quotient_digit = dividend_divisor_len_diff;

  // Initial divisor shift.
  details::FillFirstZeros(cur_divisor.digits_, dividend_length);

  // Remove signs
  result.remainder.is_negative_ = false;
  cur_divisor.is_negative_      = false;

  // Digit such that remainder - cur_digit * divisor >= 0, but remainder -
  // - (cur_digit + 1) * divisor < 0.
  DIGIT_TYPE cur_digit = 0;
  while (result.remainder >= divisor_abs_value) {
    cur_digit = FindMinimalDigitForDivision(result.remainder, cur_divisor);
    result.quotient.digits_[cur_quotient_digit--] = cur_digit;
    result.remainder -= cur_divisor * static_cast<BigInteger>(cur_digit);

    // divisor /= BASE.
    // Only if divisor is actually multiplied by BASE^n.
    if (cur_quotient_digit >= 0) {
      cur_divisor.digits_.erase(cur_divisor.digits_.begin());
    }
  }

  // Just in case.
  details::DeleteTrailingZeros(result.quotient. digits_);
  details::DeleteTrailingZeros(result.remainder.digits_);

  result.quotient.is_negative_  = is_quotient_negative;
  result.remainder.is_negative_ = is_remainder_negative; 

  return result;
}

#endif
