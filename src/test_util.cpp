// ----------------------------------------------------------------------------
//
// FILENAME: test_util.h
//
// DESCRIPTION:
//    This file used to write some test function in util to check if the
//    some functions works good
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include "poly.h"
#include "util.h"

#include <iostream>
#include <stdio.h>

int main() {
  Poly<4> p;
  printf("Test: default poly is zero: %d \n", IsZero(p));

  double coef1[4] = {3, 2, 4, 56};
  Poly<5> p1(coef1, 4);
  Monic(p1);
  std::cout << "Test: Monic polynomial: " << p1 << std::endl;

  double coef2[3] = {1, -2, 1};
  Poly<3> p2(coef2, 3);
  double coef3[2] = {-1, 1};
  Poly<3> p3(coef3, 2);
  auto gcd23(GCD(p2, p3));
  std::cout << "Test: GCD of p2 and p3: " << gcd23 << std::endl;

  std::cout << "Test: Uppper bound of p2: " << UpperBound(p2) << std::endl;

  std::cout << "Test: add 1 to p2: " << AddToX(p2, 1) << std::endl;

  double coef4[3] = {0, 1, 2};
  Poly<3> p4(coef4, 3);
  std::cout << "Test: p4 has zero root : " << ZeroRoots(&p4) << std::endl;
  std::cout << "Test: p4=" << p4 << std::endl;

  double coef5[3] = {1, -2, 1};
  Poly<3> p5(coef5, 3);
  Poly<3> sqrfree_poly5[2];
  int num_square_free5(SquareFreeDecompose(p5, sqrfree_poly5));
  printf("Test: Square free decomposition \n ");
  std::cout << "\t Original: " << p5 << std::endl;
  for (int i = 0; i < num_square_free5; i++)
    std::cout << "\t " << sqrfree_poly5[i] << std::endl;

  // Test Linear Solver
  double coef6[2] = {-445.22, 1};
  Poly<3> p6(coef6, 2);
  Range root6[1];
  int num_root_6 = 0;
  Linear(p6, 1, root6, &num_root_6);
  printf("Test: Linear Solver \n ");
  std::cout << "\t Original: " << p6 << std::endl;
  for (int i = 0; i < num_root_6; i++)
    std::cout << "\t " << root6[i].left_end << " to " << root6[i].right_end
              << std::endl;

  // Test Quadratic Solver
  double coef7[3] = {35, -12, 1};
  Poly<3> p7(coef7, 3);
  Range root7[2];
  int num_root_7 = 0;
  Quadratic(p7, 1, root7, &num_root_7);
  printf("Test: Quadratic Solver \n ");
  std::cout << "\t Original: " << p7 << std::endl;
  for (int i = 0; i < num_root_7; i++)
    std::cout << "\t " << root7[i].left_end << " to " << root7[i].right_end
              << std::endl;
  int a = 0;
  return 0;
}
