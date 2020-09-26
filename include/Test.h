#ifndef TEST_H
#define TEST_H
#include "Poly.h"

using namespace std;

class Test {
 public:
  Test(vector<Poly> polys, vector<double> ys, vector<double> gs);
  void testAccuracy();
  void testGradient();
  void testNewTown(vector<vector<double>> roots);
  void testSignRule();

  private:
  bool validRoot(double x, double target);
  vector<Poly> _polys;
  vector<double> _ys;
  vector<double> _gs;

};

#endif