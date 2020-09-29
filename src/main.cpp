#include <math.h>

#include <iostream>

#include "Poly.h"
#include "Test.h"
using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  Test tester = Test();
  tester.testAccuracy();
  tester.testGradient();
  tester.testNewton();
  tester.testSign();
  return 0;
}
