
#include <stdio.h>

#include <iostream>
#include <ostream>

#include "budan.h"
#include "poly.h"
#include "range.h"
#include "util.h"

int main() {
  // double coef[5] = {.48e-2, -.88e-1, .51, -1.2,
  // 1};  // (x - .1)(x - .3)(x - .4) ^ 2
  double coef[3] = {12, -7, 1};  // (x-3)(x-4)
  // double coef[7] = {-.8e-2, .92e-1, .674, -9.139, 12.59, -6.1, 1};
  Range roots[6];
  int num_roots = BudanRootIsolate(coef, 3, roots);

  printf("Budan's Theorem Results: %d \n", num_roots);
  for (int i = 0; i < num_roots; i++) {
    printf("left: %f, right: %f \n", roots[i].left_end, roots[i].right_end);
  }

  return 0;
}

// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 
// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 + 1e-6
// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 - 1e-6
//.3e-1*x^6-600.018*x^5+3000059.9898*x^4+1200083.99832*x^3+180015.599883*x^2+12001.139997*x+300.030000
// x^6+6.6*x^5+18.15*x^4+26.620*x^3+21.9615*x^2+9.66306*x+1.771561
// x^6+92.3521*x^2-19.22*x^4
// x^6+.3*x^5-.61*x^4-.127*x^3+.900e-1*x^2+.964e-2*x-.1680e-2
// x^6+4.2*x^4+5.88*x^2+2.744
// x^6-.196e-1+2.79*x^4+1.932*x^2
// x^6-6.1*x^5+12.59*x^4-9.139*x^3+.674*x^2+.92e-1*x-.8e-2
//- 2*x + 1 + 2*x^2
