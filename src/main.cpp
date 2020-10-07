#include <math.h>

#include <iostream>
#include <vector>

#include "coef.h"
#include "budan.h"
#include "poly.h"

using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  vector<double> coef1 = {0, 0, 1, -7.0};
  Poly test1 = Poly(coef1);
  vector<double> coef2 = {0, 0, 1, -7.001};
  Poly test2 = Poly(coef2);
  Poly test = test1* test2 * test2;

  Budan util;
  vector<Poly> ans = util.squareFreeDecompo(test);
  for(Poly tmp:ans)
  cout<< tmp<< endl;

  return 0;
}
