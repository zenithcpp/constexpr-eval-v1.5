#include "../constexpr_eval.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
  std::cout << std::setprecision(15);

  constexpr double d1 = ce::eval<"1 + 2 * 3">();
  std::cout << "1 + 2 * 3 = " << d1 << "\n";

  constexpr double d2 = ce::eval<"(1 + 2) * 3">();
  std::cout << "(1 + 2) * 3 = " << d2 << "\n";

  constexpr double d3 = ce::eval<"10 / 3">();
  std::cout << "10 / 3 = " << d3 << "\n";

  constexpr double d4 = ce::eval<"2 ^ 10">();
  std::cout << "2 ^ 10 = " << d4 << "\n";

  constexpr double d5 = ce::eval<"sin(pi/2)">();
  std::cout << "sin(pi/2) = " << d5 << "\n";

  constexpr int64_t i1 = ce::eval<"0xFF & 0b1010", ce::int_tag>();
  std::cout << "0xFF & 0b1010 = " << i1 << "\n";

  constexpr long double ld = ce::eval<"hypot(3,4)", ce::long_double_tag>();
  std::cout << "hypot(3,4) = " << static_cast<double>(ld) << "\n";

  return 0;
}
