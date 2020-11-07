#include "budan.h"
#include "vincent.h"
//#include "test.h"

using std::vector;

int main(int argc, char const *argv[]) {
  // testBudan();
  vector<double> testCoef = {1, 1, -2};
  Poly p = Poly(testCoef);
  vector<double> roots = vincentSolve(p);
  for (double root : roots) std::cout << "Root: " << root << std::endl;

  return 0;
}
