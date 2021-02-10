// ----------------------------------------------------------------------------
//
// FILENAME: test_poly.h
//
// DESCRIPTION:
//    This file used to write some test function in poly class to check if the
//    some functions works good
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include "poly.h"

#include <iostream>
#include <stdio.h>

int main() {
  // Test Poly Class

  Poly<4> p;
  std::cout << "Test: default polynomial " << p << std::endl;

  double coef1[3] = {1, -2, 1};
  Poly<3> p1(coef1, 3);
  std::cout << "Test: polynomial assign with coef " << p1 << std::endl;

  Poly<4> p2(p1);
  std::cout << "Test: polynomial assign with other poly" << p2 << std::endl;

  Poly<4> p3 = p1;
  std::cout << "Test: polynomial copy with other poly" << p3 << std::endl;

  std::cout << "Test: contain zero 4 " << p3.containZero(4) << std::endl;

  std::cout << "Test: degree of p3 " << p3.get_degree() << std::endl;

  auto v(p3.ValueAt(0));
  std::cout << "Test: p3 value at 0 from  " << v.lower() << " to " << v.upper()
            << std::endl;

  v = (p3.DerivativeAt(0));
  std::cout << "Test: p3 derivate at 0 from  " << v.lower() << " to "
            << v.upper() << std::endl;

  auto p3d(p3.Derivative());
  std::cout << "Test: derative of p3 " << p3d << " degree: " << p3d.get_degree()
            << std::endl;

  auto ld(p3.lead_coef());
  std::cout << "Test: p3 lead coef from  " << ld.lower() << " to " << ld.upper()
            << std::endl;

  std::cout << "Test: sign change of p3 " << p3.SignChange() << std::endl;

  p3 += p3;
  std::cout << "Test: p3 += p3 " << p3 << std::endl;

  p3 -= p2;
  std::cout << "Test: p3 -= p2 " << p3 << std::endl;

  p3 *= 3;
  std::cout << "Test: p3 *= 3 " << p3 << std::endl;

  p3 /= 3;
  std::cout << "Test: p3 /= 3 " << p3 << std::endl;

  auto p4 = p3 + p3;
  std::cout << "Test: p3 + p3 " << p4 << std::endl;

  auto p5 = p4 - p3;
  std::cout << "Test: p4 - p3 " << p5 << std::endl;

  auto p6 = p3 * p3;
  std::cout << "Test: p3 * p3 " << p6 << std::endl;

  double coef7[2] = {-1, 1};
  Poly<2> tmp(coef7, 2);
  auto p7 = Division(p2, tmp);
  std::cout << "Test: p2 * tmp: quotient " << p7.quotient
            << " remainder=" << p7.remainder << std::endl;
  return 0;
}
