#include <stdio.h>

#include "poly.h"

int main() {
  double coef1[3] = {2.0, -3.0, 1.0};  // x^2 - 2x + 1
  Poly<3> p1(coef1, 3);

  double coef2[2] = {-1.0, 1.0};  // x - 1
  Poly<3> p2(coef2, 2);

  auto ans(p1 / p2);

  std::cout << "Division ans : \n "
            << "Qutient: " << ans.quotient << "\n Remainder : " << ans.remainder
            << std::endl;

  return 0;
}
