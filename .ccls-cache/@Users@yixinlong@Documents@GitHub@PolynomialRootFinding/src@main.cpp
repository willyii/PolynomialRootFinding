#include "budan.h"
#include "vincent.h"
//#include "test.h"

using std::vector;

int main(int argc, char const *argv[]) {
  // testBudan();
  vector<double> testCoef = {1, -0.8, 0.15};
  Poly p = Poly(testCoef);
  vector<double> roots = vincentSolve(p);
  for (double root : roots) std::cout << "Testd: " << root << std::endl;
  return 0;
}
//
//

//#include <gsl/gsl_poly.h>
//#include <stdio.h>

// int main(void) {
//  int i;
//  [> coefficients of P(x) =  -1 + x^5  <]
//  double a[6] = {-1, 0, 0, 0, 0, 1};
//  double z[10];

//  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(6);

//  gsl_poly_complex_solve(a, 6, w, z);

//  gsl_poly_complex_workspace_free(w);

//  for (i = 0; i < 5; i++) {
//    printf("z%d = %+.18f %+.18f\n", i, z[2 * i], z[2 * i + 1]);
//  }

//  return 0;
//}
