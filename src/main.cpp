#include <math.h>

#include <iostream>

#include "Poly.h"
#include "Test.h"
using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  cout << "Hello World" << endl;
  cout << "Fuck you" << endl;
  Poly tmpPoly;
  vector<Poly> testPolys;
  vector<double> ys, gs;
  vector<vector<double>> roots;

  // /* x^2 + 2x + 1 || f(3) = 16 || no root */
  // tmpPoly = Poly(vector<double>{1.0, 2.0, 1.0});
  // testPolys.push_back(tmpPoly);
  // ys.push_back(16);
  // gs.push_back(8);
  // roots.push_back(vector<double>{});

  /* x^2 - 8x + 15 || f(3) = 0 || root in {3, 5} */
  tmpPoly = Poly(vector<double>({1.0, -8.0, 15.0}));
  testPolys.push_back(tmpPoly);
  ys.push_back(0);
  gs.push_back(-2);
  roots.push_back(vector<double>{3, 5});

  /* x^3 - 12x^2 + 47x - 60 || f(3) = 36 || root in {3,4,5}*/
  tmpPoly = Poly(vector<double>({1, -12.0, 47.0, -60.0}));
  testPolys.push_back(tmpPoly);
  ys.push_back(0);
  gs.push_back(2);
  roots.push_back(vector<double>{3, 4, 5});

  Test tester(testPolys, ys, gs);
  tester.testAccuracy();
  tester.testGradient();
  tester.testNewTown(roots);
  tester.testSignRule();
  return 0;
}
