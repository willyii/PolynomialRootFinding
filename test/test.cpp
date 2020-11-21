#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "budan.h"
#include "gsl/gsl_poly.h"
#include "vincent.h"

#define PRINTROOT true

bool validAnswer(vector<double> actual, vector<double> current) {
  bool state = true;
  if (actual.size() != current.size()) state = false;
  if (state == true) {
    for (size_t i = 0; i < actual.size(); i++)
      if (fabs(actual[i] - current[i]) > TESTERROR) state = false;
  }

  if (PRINTROOT) {
    std::cout << "Acutal roots: ";
    for (size_t i = 0; i < actual.size(); i++) std::cout << " " << actual[i];
    std::cout << "\n";
    std::cout << "Current roots: ";
    for (size_t i = 0; i < current.size(); i++) std::cout << " " << current[i];
    std::cout << "\n";
  }
  return state;
}

vector<double> gslSolve(double coef[], int N) {
  // ------------- State of Art Result -------------------
  double gslAnsComplex[(N - 1) * 2];
  vector<double> gsl_ans;
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(N);
  gsl_poly_complex_solve(coef, N, w, gslAnsComplex);
  gsl_poly_complex_workspace_free(w);

  for (int i = 0; i < (N - 1) * 2; i++) {
    if (gslAnsComplex[2 * i + 1] == 0 &&
        fabs(gslAnsComplex[2 * i] - 0) > EPSILON) {
      gsl_ans.emplace_back(gslAnsComplex[2 * i]);
    }
  }
  sort(gsl_ans.begin(), gsl_ans.end());
  return gsl_ans;
}

int rand_int(int b = 5) { return rand() % b + 2; }

double rand_float(double a = -10, double b = 10) {
  return ((double)rand() / RAND_MAX) * (b - a) + a;
}

/* Process polynomial from command line */
bool validSinglePoly(double coef[], int N) {
  bool state = true;
  vector<double> gsl_ans = gslSolve(coef, N);

  std::reverse(coef, coef + N);
  vector<double> coefv;
  coefv.assign(coef, coef + N);
  Poly testPoly(coefv);
  std::cout << testPoly << std::endl;

  std::cout << "----"
            << "Validating budan  "
            << "-----" << std::endl;
  vector<double> budan_ans = budanSolve(testPoly);
  std::sort(budan_ans.begin(), budan_ans.end());
  state = validAnswer(gsl_ans, budan_ans) && state;

  std::cout << "----"
            << "Validating vincent"
            << "-----" << std::endl;
  vector<double> vincent_ans = vincentSolve(testPoly);
  std::sort(vincent_ans.begin(), vincent_ans.end());
  state = validAnswer(gsl_ans, vincent_ans) && state;

  return state;
}

/* Porcess random generated polynomial */
int testRandomPoly() {
  std::cout << "===============================" << std::endl;
  std::cout << "Test a single random polynomial:" << std::endl;
  // int N = rand_int();
  // double coef[N];
  // for (size_t i = 0; i < N; i++) {
  //  coef[i] = rand_float();
  //}

  /* TODO: Vincent zero root */
  // x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 
  // x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 + 1e-6
  // x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 - 1e-6
  // .3e-1*x^6-600.018*x^5+3000059.9898*x^4+1200083.99832*x^3+180015.599883*x^2+12001.139997*x+300.030000
  // x^6+6.6*x^5+18.15*x^4+26.620*x^3+21.9615*x^2+9.66306*x+1.771561
  // x^6+92.3521*x^2-19.22*x^4
  // x^6+.3*x^5-.61*x^4-.127*x^3+.900e-1*x^2+.964e-2*x-.1680e-2
  // x^6+4.2*x^4+5.88*x^2+2.744
  // x^6-.196e-1+2.79*x^4+1.932*x^2
  // x^6-6.1*x^5+12.59*x^4-9.139*x^3+.674*x^2+.92e-1*x-.8e-2

  int N = 7;
  double coef[] = {1, -6.1, +12.59, -9.139, +.674, +.92e-1, -.8e-2};
  // double coef[] = {1, 0, 2.79, 0, 1.932, 0, -.196e-1};
  // double coef[] = {1, 0, +4.2, 0, +5.88, 0, +2.744};
  // double coef[] = {1, +.3, -.61, -.127, +.900e-1, +.964e-2, -.1680e-2};
  // double coef[] = {1, 0, -19.22, 0, +92.3521, 0, 0};
  // double coef[] = {.3e-1,         -600.018,
  //                 +3000059.9898, +1200083.99832 + 180015.599883,
  //                 +12001.139997, +300.030000};
  // double coef[] = {1, +6.6, +18.15, +26.620, +21.9615, 9.66306,
  // +1.771561};
  std::reverse(coef, coef + N);

  bool result = validSinglePoly(coef, N);
  if (result) {
    std::cout << "Test on random polynomial passed\n" << std::endl;
  } else {
    std::cout << "Test on random polynomial failed\n" << std::endl;
  }

  return 0;
}

/* TODO: Test Polynomials from file */
int testFromFile(char const *filename) { return 0; }

int testSinglePoly(int argc, char const *argv[]) {
  std::cout << "===============================" << std::endl;
  std::cout << "Test a polynomial from input  :" << std::endl;
  int N = argc - 1;
  double coef[N];
  for (size_t i = 0; i < N; i++) coef[i] = std::stod(argv[i + 1]);
  std::reverse(coef, coef + N);
  bool result = validSinglePoly(coef, N);
  if (result) {
    std::cout << "Test on input polynomial passed\n" << std::endl;
  } else {
    std::cout << "Test on input polynomial failed\n" << std::endl;
  }

  return 0;
}

int main(int argc, char const *argv[]) {
  int state;
  srand(time(NULL));
  if (argc == 1) {
    state = testRandomPoly();
    if (state == 1) {
      std::cout << "Test random polynomia function error\n" << std::endl;
    }
    return state;
  } else if (argc == 2) {
    state = testFromFile(argv[1]);
    if (state == 1) {
      std::cout << "Test from file function error\n" << std::endl;
    }
    return state;
  } else {
    state = testSinglePoly(argc, argv);
    if (state == 1) {
      std::cout << "Test single polynomial function error\n" << std::endl;
    }
    return state;
  }
  return 0;
}
