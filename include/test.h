#ifndef TEST_H
#define TEST_H
#include <stdlib.h>
#include <time.h>

#include <cassert>

#include "budan.h"

using namespace std;

// random float generator
double rand_float(double a = -10, double b = 10) {
  return ((double)rand() / RAND_MAX) * (b - a) + a;
}

void testPoly(Poly& p, Budan util) {
  vector<double> roots;
  // cout << "===============================" << endl;
  cout << "Poly: " << p << endl;
  roots = util.solve(p);
  cout<<"There are "<<roots.size()<< " roots in total" << endl;
  for (auto root : roots) {
    cout << "root: " << root <<" Value: "<<p.valueAt(root) << endl;
    assert(p.valueAt(root) == 0.0);
  }
  cout << "Test Pass" << endl;
  cout << "===============================" << endl;
}

void testBudan() {
  srand(time(NULL));
  vector<double> tmpCoef, roots;
  Poly test1, test2, test3, test;
  Budan util;

  // 1 root
  // tmpCoef = {1, -2, 1};
  // test = Poly(tmpCoef);
  // testPoly(test, util);

  // repeat roots
  tmpCoef = {1, rand_float()};
  test1 = Poly(tmpCoef);
  tmpCoef = {1, rand_float()};
  test2 = Poly(tmpCoef);

  tmpCoef = {1, rand_float()};
  test3 = Poly(tmpCoef);
  test = test1 * test2 *test3;
  testPoly(test, util);

  // Random poly
  // tmpCoef = {3.02385, 1.91351, 0.288506, 8.9231};
  // tmpCoef = { rand_float(), rand_float(), rand_float(), rand_float(),  rand_float()};
  // test = Poly(tmpCoef);
  // // cout<<test.valueAt(1.50220411)<<endl;
  // testPoly(test, util);
}

#endif