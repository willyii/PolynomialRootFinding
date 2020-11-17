#include <stdlib.h>
#include <time.h>

#include "budan.h"
#include "gsl/gsl_poly.h"
#include "vincent.h"

using std::cout;
using std::endl;
using std::reverse;
using std::sort;

// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 
// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 + 1e-6
// x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 - 1e-6

// random float generator
double rand_float(double a = -10, double b = 10) {
  return ((double)rand() / RAND_MAX) * (b - a) + a;
}

// random float generator
int rand_int(int b = 5) { return rand() % b + 2; }

vector<double> gslSolve(double coef[], int N) {
  // ------------- State of Art Result -------------------
  double gslAnsComplex[(N - 1) * 2];
  vector<double> gsl_ans;
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(N);
  gsl_poly_complex_solve(coef, N, w, gslAnsComplex);
  gsl_poly_complex_workspace_free(w);

  for (int i = 0; i < (N - 1) * 2; i++) {
    if (gslAnsComplex[2 * i + 1] == 0 && gslAnsComplex[2 * i] != 0) {
      gsl_ans.emplace_back(gslAnsComplex[2 * i]);
    }
  }
  sort(gsl_ans.begin(), gsl_ans.end());
  return gsl_ans;
}

int main(int argc, char const *argv[]) {
  srand(time(NULL));

  // ------------- Random Poly     -------------------
  // int N = rand_int();
  //// int N = 6;
  // double coef[N];
  // for (size_t i = 0; i < N; i++) {
  //  coef[i] = rand_float();
  //}
  double coef[] = {1, -1.2, .51, -.88e-1, .48e-2 - 1e-6};
  std::reverse(coef, coef + 5);
  int N = 5;
  vector<double> gsl_ans = gslSolve(coef, 5);

  vector<double> coefVec;
  coefVec.assign(coef, coef + N);
  reverse(coefVec.begin(), coefVec.end());
  Poly testPoly = Poly(coefVec);
  // ------------- Budan Result     -------------------
  vector<double> budan_ans = budanSolve(testPoly);
  sort(budan_ans.begin(), budan_ans.end());

  // ------------- Budan Result     -------------------
  vector<double> vin_ans = vincentSolve(testPoly);
  sort(vin_ans.begin(), vin_ans.end());

  // check the root size
  cout << "Test Polynomial: " << testPoly << endl;
  cout << "=====================================" << endl;
  cout << "Vlaidating the correctness # of roots" << endl;
  if (budan_ans.size() != gsl_ans.size())
    cout << "Budan Test Failed, actual root size: " << gsl_ans.size()
         << ", Budan root size: " << budan_ans.size() << endl;
  if (vin_ans.size() != gsl_ans.size()) {
    cout << "Vincent Test Failed, actual root size: " << gsl_ans.size()
         << ", Vincent root size: " << vin_ans.size() << endl;
  }
  cout << "Vlaidation of # of roots finished\n" << endl;

  // Checkt the coorecness of rot
  cout << "=====================================" << endl;
  cout << "Vlaidating the correctness of roots" << endl;
  cout << "Acutal Roots: " << endl;
  for (int i = 0; i < gsl_ans.size(); i++) {
    cout << gsl_ans[i] << ", ";
  }
  cout << "\n";
  for (int i = 0; i < budan_ans.size(); i++) {
    if (fabs(gsl_ans[i] - budan_ans[i]) > TESTERROR)
      cout << "Budan test failed in one root: actual value " << gsl_ans[i]
           << " , my version: " << budan_ans[i] << endl;
  }
  for (int i = 0; i < vin_ans.size(); i++) {
    if (fabs(gsl_ans[i] - vin_ans[i]) > TESTERROR)
      cout << "Vincent test failed in one root: actual value " << gsl_ans[i]
           << " , my version: " << vin_ans[i] << endl;
  }
  cout << "Vlaidation of correctness finished\n" << endl;

  return 0;
}
