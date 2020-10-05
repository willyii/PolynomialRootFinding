#include <math.h>

#include <iostream>
#include <vector>

#include "poly.h"

using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  // Test tester = Test();
  // tester.testAccuracy();
  // tester.testGradient();
  // tester.testNewton();
  // tester.testSign();
  // tester.testGcd();
  vector<double> x = {1, 0, -1};
  Poly tmp1 = Poly(x);
  vector<double> y = {0, 1, -1};
  Poly tmp2 = Poly(y);
  Poly tmp3 = tmp1 / tmp2;
  cout << "Tmp1: " << tmp1.getGradPoly() << endl;
  cout << "Value at 1: "<< tmp1.valueAt(1)<<endl;
  cout << "Value at -1: "<< tmp1.valueAt(-1)<<endl;


  return 0;
}
