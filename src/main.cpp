#include <math.h>

#include <iostream>
#include <vector>

#include "coef.h"
#include "budan.h"
#include "poly.h"

using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  // vector<double> coe = {1, -4, -122, 252, 3969};
  vector<double> coe = {1, 2, 1};
  Budan util;
  Poly test = Poly(coe), tmp = util.addToP(test, -2);
  int ans = util.signChangeNum( tmp );
  cout << "test: " << ans << endl;

  return 0;
}
