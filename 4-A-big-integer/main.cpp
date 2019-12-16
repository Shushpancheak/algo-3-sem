#include "biginteger.h"
#include <iostream>


int main() {

  const size_t len = 100;

  BigInteger arr[len];

  for (size_t i = 0; i < len; ++i) {
    arr[i] = (i + 123) * 8135 * pow(-1, i);
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;

  for (size_t i = 1; i < len; ++i) {
    arr[i] = arr[i - 1] * arr[i] * 1234 * (-1);
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;

  for (size_t i = 1; i < len; ++i) {
    arr[i] = arr[i] / arr[i - 1];
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;

  for (size_t i = 1; i < len; ++i) {
    arr[i] = arr[i - 1] * arr[i] * 1234 * (-1);
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;


  for (size_t i = 1; i < len; ++i) {
    arr[i] = arr[i - 1] % arr[i];
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;


  for (size_t i = 1; i < len; ++i) {
    arr[i] = arr[i - 1] + arr[i - 1];
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;

  return 0;

  BigInteger a1, a2;
  std::cin >> a1 >> a2;

  std::cout << a2 * a1;


  return 0;

  BigInteger big = 234567;

  big *= big;
  big *= big;
  big *= big;
  big *= big;
  big *= big;
  big *= big;
  big *= big;
  big *= big;
  big *= big;
  std::cout << big * big;

  //return 0;


  BigInteger b;
  std::cin >> b;
  b *= b;
  b /= b;
  b += b;
  b -= b;


  int64_t number = 999;
  int iterations = 10;

  BigInteger a = 10;
  std::cout << a << std::endl;

  for (int i = -iterations; i < iterations; ++i) {
    a += i;
    a *= i;
    a += a * a;
    std::cout << a << std::endl;
  }

  std::cout << "new value = ";
  //std::cin >> a;
  std::cout << "So, new value is: " << a << std::endl;
  for (size_t i = 0; i < iterations; ++i) {
    BigInteger quo = a / (number + 1);
    BigInteger rem = a % number;

    std::cout << i << " quo: " << quo << std::endl;
    std::cout << i << " rem: " << rem << std::endl;

    a = quo;
  }

  std::cout << "==== " << BigInteger(9890) / 1000 << std::endl;
  std::cout << "=====" << a - 123089721 << " " << a++ << " " << ++a << " " << a-- << " " << --a << " " << a.toString() << " " << bool(a) << " " << int(a) << " " << a / 1234125 << " " << a * 12341235 << " " << a +  21234234 << " " << a * 123412352 << " " << a % 3123;

  std::cout << " | " << BigInteger(123456789) * 123456789 * 123456789;

  return 0;
}