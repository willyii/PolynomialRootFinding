#include <iostream>

#include "poly.h"

int main() {
  // Poly a = Poly();
  double coef[3] = {1.0, -2.0, 1.0};
  Poly<3> a = Poly<3>(coef, 3);

  std::cout << a.gradientAt(1) << std::endl;
  return 0;
}
