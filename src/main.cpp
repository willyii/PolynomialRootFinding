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

  interval a = {-1.0, 2.0};
  interval b = a;
  printf("= operator b left %f b right %f \n", b.left, b.right);
  b += a;
  printf("+ operator b left %f b right %f \n", b.left, b.right);
  b -= a;
  printf("- operator b left %f b right %f \n", b.left, b.right);
  b *= a;
  printf("* operator b left %f b right %f \n", b.left, b.right);
  b /= a;
  printf("/ operator b left %f b right %f \n", b.left, b.right);

  printf("====================================================\n");
  interval c = a + a;
  printf("+ operator c left %f b right %f \n", c.left, c.right);
  c *= -9;
  printf("* operator c left %f b right %f \n", c.left, c.right);
  c = 3 - c;
  printf("- operator c left %f b right %f \n", c.left, c.right);
  c = c / -3;
  printf("/ operator c left %f b right %f \n", c.left, c.right);
  c = 3 / c;
  printf("/ operator c left %f b right %f \n", c.left, c.right);

  double coef[3] = {1, -2, 1};
  Poly<4> p1(coef, 3);
  double coef2[2] = {-1, 1};
  Poly<3> p2(coef, 2);
  std::cout << Quotient(p1, p2) << std::endl;

  return 0;
}
