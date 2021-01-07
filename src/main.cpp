#include <iostream>

//#include "poly.h"
#include "util.h"

int main() {
  // Poly a = Poly();
  double coef[2] = {-9.0, 1.0};  // x^2 - 2 x + 1
  Poly<2> a(coef, 2);
  double coef2[3] = {1.0, -2.0, 1.0};  // x^2 - 2 x + 1
  Poly<2> b(coef2, 3);
  // Poly<1> a(b);
  auto c = b / a;

  std::cout << c.q << std::endl;
  std::cout << c.r << std::endl;
  std::cout << upperBound(b) << std::endl;
  std::cout << lowerBound(b) << std::endl;
  std::cout << signChangeNum(b) << std::endl;
  std::cout << addToX(b, 2) << std::endl;
  // printf("N of c: %d\n", c.size());
  return 0;
}
