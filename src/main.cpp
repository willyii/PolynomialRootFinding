#include <math.h>

#include <iostream>
#include <vector>
#include <time.h>

#include "coef.h"
#include "budan.h"
#include "poly.h"

using namespace std;


double rand_float(double a= -100, double b=100) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}


int main(int argc, char const *argv[]) {
  srand(time(NULL)) ;
  /* code */
  vector<double> coef1 = {0, 0, 1, -2, -23, 60};
  Poly test1 = Poly(coef1);

  Budan util;
  vector<double> roots = util.solve(test1);
  for(auto ans:roots)
  cout<<"root: "<<ans<< "\n";


  return 0;
}
