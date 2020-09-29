#include "Test.h"

#include <math.h>

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "Budan.h"
#include "Param.h"
#include "Poly.h"

using namespace std;

/* Test the accuracy of the poly*/
void Test::testAccuracy() {
  cout << "Testing Accuracy... " << endl;
  Poly testPoly;

  // Case 1, x^2 + 2x + 1
  testPoly = Poly(vector<double>{1, 2, 1});
  assert(testPoly.getValue(1) == 4);
  assert(testPoly.getValue(2) == 9);

  // Case 2, x^3 - 11x^2 + 39x -45
  testPoly = Poly(vector<double>{1, -11, 39, -45});
  assert(testPoly.getValue(3) == 0);
  assert(testPoly.getValue(5) == 0);

  cout << "Test Accuracy Success" << endl;
}

/* Test the gradient function*/
void Test::testGradient() {
  cout << "Testing Gradient... " << endl;
  Poly testPoly;

  // Case 1, x^2 + 2x + 1 || g(1) = 4, g(2) = 6
  testPoly = Poly(vector<double>{1, 2, 1});
  assert(testPoly.getGradient(1) == 4);
  assert(testPoly.getGradient(2) == 6);

  // Case 2, x^3 - 11x^2 + 39x -45 || g(3) = 0, g(5) = 4
  testPoly = Poly(vector<double>{1, -11, 39, -45});
  assert(testPoly.getGradient(3) == 0);
  assert(testPoly.getGradient(5) == 4);

  cout << "Test Gradient Success " << endl;
}

/* Test Newton Ralphson function*/
void Test::testNewton() {
  cout << "Testing Newton... " << endl;
  Poly testPoly;
  Budan testBudan = Budan();
  double testRoot;

  // Case1, x^2 - 2x + 1 || root = 1
  testPoly = Poly(vector<double>{1, -2, 1});
  testRoot = testBudan.NewtonRaphson(testPoly, 3.0);
  assert(abs(testPoly.getValue(testRoot) - 0.0) <= Param::EPSILON);

  // Case2, x^3 - 11x^2 + 39x -45 || root =3 or root = 5
  testPoly = Poly(vector<double>{1, -11, 39, -45});
  testRoot = testBudan.NewtonRaphson(testPoly, 7.0);
  assert(abs(testPoly.getValue(testRoot) - 0.0) <= Param::EPSILON);
  testRoot = testBudan.NewtonRaphson(testPoly, 1.0);
  assert(abs(testPoly.getValue(testRoot) - 0.0) <= Param::EPSILON);

  cout << "Test Newton Success " << endl;
}

/* Test Sign change rule*/
void Test::testSign() {
  cout << "Testing SignRule... " << endl;
  Poly testPoly;
  Budan testBudan = Budan();
  int left, right;

  // Case 1 (x-3)(x-4)(x-5)
  testPoly = Poly(vector<double>{1, -12, 47, -60});

  left = testBudan.signChangeNums(testPoly, 1);
  right = testBudan.signChangeNums(testPoly, 2);
  assert(left - right == 0);

  left = testBudan.signChangeNums(testPoly, 2.5);
  right = testBudan.signChangeNums(testPoly, 3.5);
  assert(left - right == 1);

  left = testBudan.signChangeNums(testPoly, 3.5);
  right = testBudan.signChangeNums(testPoly, 4.5);
  assert(left - right == 1);

  // Case 2 (x-3)^2(x-5)
  testPoly = Poly(vector<double>{1, -11, 39, -45});

  left = testBudan.signChangeNums(testPoly, 2.5);
  right = testBudan.signChangeNums(testPoly, 3.5);
  assert(left - right == 2);

  left = testBudan.signChangeNums(testPoly, 4.5);
  right = testBudan.signChangeNums(testPoly, 5.5);
  assert(left - right == 1);

  cout << "Testing SignRule Pass " << endl;
}
