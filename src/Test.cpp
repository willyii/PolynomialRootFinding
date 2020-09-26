#include "Test.h"

#include <math.h>

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

Test::Test(vector<Poly> polys, vector<double> ys, vector<double> gs)
    : _polys(polys), _ys(ys), _gs(gs) {}

/* Test if this poly can get correct value*/
void Test::testAccuracy() {
  cout << "=================================" << endl;
  cout << "========Test the correctness=====" << endl;
  int N = _polys.size();
  for (int i = 0; i < N; i++) {
    if (_polys[i].getValue(3.0) != _ys[i]) {
      cout << "==============Test fialed========" << endl;
      cout << "Expected Value: " << _ys[i] << endl;
      cout << "Actual Value: " << _polys[i].getValue(3.0) << endl;
      return;
    }
  }
  cout << "===========Test sucess==========="
       << "\n\n"
       << endl;
  return;
}

/* Test if this poly can get correct gradient*/
void Test::testGradient() {
  cout << "=================================" << endl;
  cout << "========Test the gradient=====" << endl;
  int N = _polys.size();
  for (int i = 0; i < N; i++) {
    if (_polys[i].getGradient(3.0) != _gs[i]) {
      cout << "==============Test fialed========" << endl;
      cout << "Expected Value: " << _gs[i] << endl;
      cout << "Actual Value: " << _polys[i].getValue(3.0) << endl;
      return;
    }
  }
  cout << "===========Test sucess==========="
       << "\n\n"
       << endl;
  return;
}

/* Test if NewtownRaphs can work correctly*/
void Test::testNewTown(vector<vector<double>> roots) {
  cout << "=================================" << endl;
  cout << "========Test the NewTown=====" << endl;

  double epsilon = numeric_limits<double>::epsilon();
  int N = _polys.size();
  double testroot;
  bool pass;

  for (int i = 0; i < N; i++) {
    testroot = _polys[i].NewtonRaphson(7.0);
    pass = false;
    for (auto r : roots[i]) {
      pass = pass || validRoot(testroot, r);
    }
    if (!pass) {
      cout << "==============Test fialed========" << endl;
      cout << "Case " << i << "  Test Root: " << testroot << endl;
      return;
    }
  }
  cout << "===========Test sucess==========="
       << "\n\n"
       << endl;
  return;
}

bool Test::validRoot(double x, double y) {
  return abs(x - y) <= Param::EPSILON;
}

/* Test the correctness of Sign Rule*/
void Test::testSignRule() {
  cout << "==========Test Sign Change Count=========" << endl;
  vector<double> test1;
  Poly testPoly;
  int changes;

  // Case1
  test1 = {1, -12.0, 47.0, -60.0};
  testPoly = Poly(test1);
  changes = testPoly.signChangeNums(1);
  assert(changes == 3);

  // Case2
  test1 = {1, -12.0, 47.0, -60.0};
  testPoly = Poly(test1);
  changes = testPoly.signChangeNums(-4);
  assert(changes == 3);

  // Case3
  test1 = {3, 0, -7.0, -7.0};
  testPoly = Poly(test1);
  changes = testPoly.signChangeNums(0) - testPoly.signChangeNums(2);
  assert(changes == 2);

  // Case4
  test1 = { 1, -2.0, 1.0};
  testPoly = Poly(test1);
  changes = testPoly.signChangeNums(0) - testPoly.signChangeNums(8);
  cout<<"Changes : "<<changes<<endl;
  // assert(changes == 2);

  cout << "==========Test Sign Change Pass=========" << endl;
}