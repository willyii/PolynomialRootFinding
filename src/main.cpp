// ----------------------------------------------------------------------------
//
// FILENAME: main.cpp
//
// DESCRIPTION:
//    This file define the program entrance of this project. It takes several
//    parameters
//
// PARAMETERS:
//    No parameter: test random generated polynomials, print running time
//    One parameter: Usually take polynomial file path as parameter. If use
//                   "valid" as  parameter, it will test random polynomials and
//                   print out the results
//    Two parameters:First parameter should be the path to poly file.
//                   If second parameter is "valid", it will print out the
//                   result get from two methods.
//                   Else second parameter will be regard as path to solution of
//                   first file
//
//
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------
#include "interval.h"
#include "poly.h"
#include "util.h"
//#include "param.h"
//#include "testfunction.h"

#include <iomanip>
#include <ios>
#include <string>

/**
 * Return true if file exists
 */
bool is_file_exist(const char *fileName) {
  std::ifstream infile(fileName);
  return infile.good();
}

int main(int argc, char *argv[]) {
  std::cout << std::setprecision(16);
  // if (argc == 1) {
  /**
   * Generate random polynomials and write to file
   * Then read from file and record the running time
   */
  //  RandomPolyToFile(kNUMPOLY);
  //  RunningTimeTest(kRANDOM_FILE);
  //} else if (argc == 2) {
  //  if (std::string(argv[1]) == "valid") {
  /**
   * Generate random polynomials and corresponding results
   * Save to file then read and slove, print out the result
   */
  //    RandomPolyToFile(kNUMPOLY);
  //    RunningCorrectTest(kRANDOM_FILE, kRANDOM_FILE_SOL);
  //    return 0;
  //  }

  /**
   * Test running time of polynomals in files
   */
  //  RunningTimeTest(argv[1]);
  //} else if (argc == 3) {
  /**
   * Test with specific file, and print out the result
   */
  //  if (std::string(argv[2]) == "valid") {
  //    RunningCorrectTest(argv[1]);
  //    return 0;
  //  }

  /**
   * Test with specific file and corresponding solution file
   */
  //  RunningCorrectTest(argv[1], argv[2]);

  //  return 0;
  //} else {
  //  printf("Invalid parameters\n");
  //  return 1;
  //}

  double coef[3] = {1, -2, 1};
  Poly<4> p1(coef, 3);
  double coef2[2] = {-1, 1};
  Poly<3> p2(coef2, 2);
  std::cout << "Test: P2 = " << p2 << std::endl;

  Poly<5> p3(p2);
  std::cout << "Test: P3 = " << p3 << std::endl;

  Poly<5> p4;
  p4 = p2;
  std::cout << "Test: P4 = " << p4 << std::endl;

  std::cout << "Test: P4(x=2) = " << p4.ValueAt(2) << std::endl;

  std::cout << "Test: P4 derivate at 2 = " << p4.DerivativeAt(2) << std::endl;

  std::cout << "Test: P4 derivate = " << p4.Derivative() << std::endl;

  p4 += p1;
  std::cout << "Test: P4 + p1 " << p4 << std::endl;

  p4 += 1;
  std::cout << "Test: P4 + 1 " << p4 << std::endl;

  p4 -= 1;
  std::cout << "Test: P4 - 1 " << p4 << std::endl;

  p4 -= p1;
  std::cout << "Test: P4 - p1 " << p4 << std::endl;

  p4 /= 2;
  std::cout << "Test: P4 / 2  " << p4 << std::endl;

  p4 *= 2;
  std::cout << "Test: P4 * 2  " << p4 << std::endl;

  auto p5 = p4 + p1;
  std::cout << "Test: P5  " << p5 << std::endl;

  auto p6 = p5 - p1;
  std::cout << "Test: P6  " << p6 << std::endl;

  auto p7 = p4 * p4;
  std::cout << "Test: P7  " << p7 << std::endl;

  auto p8 = p4 * 4;
  std::cout << "Test: P8  " << p8 << std::endl;

  auto p9 = Division(p1, p2);
  std::cout << "Test: p1/p2 \n "
            << " quotient: " << p9.quotient << " \n remainder: " << std::endl;

  auto p10 = Quotient(p1, p2);
  std::cout << "Test: p10 " << p10 << std::endl;

  auto p11 = Remainder(p1, p2);
  std::cout << "Test: p11 " << p11 << "  IsZero: " << IsZero(p11) << std::endl;

  auto p12 = GCD(p1, p2);
  std::cout << "Test: p12 " << p12 << std::endl;

  double test_coef[kMAXDEGREE] = {.48e-2 + 1e-6, -.88e-1, .51, -1.2, 1};
  Poly<kMAXDEGREE> poly(test_coef, 5);
  Poly<kMAXDEGREE> square_free_polys[kMAXDEGREE];
  int num_square = SquareFreeDecompose(poly, square_free_polys);

  for (int i = 0; i < num_square; i++) {
    std::cout << "Test: Square Free Polynomial: " << square_free_polys[i]
              << std::endl;
  }

  std::cout << "Test: upperbound of poly " << UpperBound(poly) << std::endl;

  std::cout << "Test: add 1 to x " << AddToX(poly, 1) << std::endl;

  return 0;
}
