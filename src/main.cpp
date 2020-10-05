#include <math.h>

#include <iostream>
#include <vector>

#include "coef.h"
#include "poly.h"

using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  Poly tmp1, tmp2, tmp3;
  vector<double> x = {1, 0, -1};
  tmp1 = Poly(x);
  vector<double> y = {0, 1, -1};
  tmp2 = Poly(y);
  tmp3 = tmp1 / tmp2;
  cout << "Tmp1: " << tmp1.getGradPoly() << endl;
  cout << "Value at 1: "<< tmp1.valueAt(1)<<endl;
  cout << "Value at -1: "<< tmp1.valueAt(-1)<<endl;


  return 0;
}
