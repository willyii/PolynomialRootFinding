#include <stdio.h>

#include "poly.h"
#include "util.h"

int main() {
  double coef1[3] = {2.0, -3.0, 1.0};  // x^2 - 2x + 1
  Poly<3> p1(coef1, 3);

  double coef2[2] = {-1.0, 1.0};  // x - 1
  Poly<3> p2(coef2, 2);

  auto ans = GCD(p2, p1);

  std::cout << "GCD ans : " << ans << std::endl;

  return 0;
}
