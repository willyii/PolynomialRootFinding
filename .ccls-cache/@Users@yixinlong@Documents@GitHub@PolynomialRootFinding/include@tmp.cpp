#include <stdio.h>

#include <iostream>

//#include "poly.h"
#include "util.h"

int main() {
  // Poly a = Poly();
  double coef[2] = {-9.0, 1.0};  // x^2 - 2 x + 1
  Poly<2> a = Poly<2>(coef, 2);
  double coef2[3] = {1.0, -2.0, 1.0};  // x^2 - 2 x + 1
  Poly<2> b = Poly<2>(coef2, 3);
  auto c = b / a;

  std::cout << c.q << std::endl;
  std::cout << c.r << std::endl;
  std::cout << upperBound(b) << std::endl;
  // printf("N of c: %d\n", c.size());
  return 0;
}
