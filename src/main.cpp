#include "bisect.h"
#include "budan.h"
#include "test.h"

using std::vector;

int main(int argc, char const *argv[]) {
  /* code */
  // testBudan();
  vector<double> testCoef = {1, -8, 15};
  Poly p = Poly(testCoef);
  vector<double> roots = bisectSolve(p);
  for (double root : roots) std::cout << "Testd: " << root << std::endl;
  return 0;
}
