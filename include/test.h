#ifndef TEST_H
#define TEST_H
#include <stdlib.h>
#include <time.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include "budan.h"
#include "vincent.h"

using namespace std;

// random float generator
double rand_float(double a = -10, double b = 10) {
  return ((double)rand() / RAND_MAX) * (b - a) + a;
}

bool validSinglePoly(Poly& p, Vincent& util, vector<double>& ans) {
  vector<double> roots = util.solve(p);
  if (ans.size() != roots.size()) {
    cout << "Validation fail on: " << p << "\n"
         << "Current root num: " << roots.size() << "\n"
         << "Actual root num: " << ans.size() << "\n";

    cout << "Current roots : "
         << "\n";
    for (auto r : roots) cout << r << "\t";
    cout << "\n";

    cout << "Actual roots : "
         << "\n";
    for (auto r : ans) cout << r << "\t";
    cout << "\n";
    return false;
  }
  sort(roots.begin(), roots.end());

  for (int i = 0; i < ans.size(); i++) {
    if (fabs(roots[i] - ans[i]) > TESTERROR) {
      cout << "Validation fail on " << p << "\n"
           << " with test root: " << roots[i] << " and actual root: " << ans[i]
           << "\n"
           << " Error: " << fabs(roots[i] - ans[i]) << endl;
      return false;
    }
  }
  return true;
}

void validPolyFromFile(string path) {
  Vincent util;
  ifstream validfile(path);
  string line, tmp_ans, tmp_ceof;
  vector<double> coef, ans;
  Poly testPoly;
  int test_count = 0, pass = 0;

  while (getline(validfile, line)) {
    coef = {};
    ans = {};
    istringstream iss(line);
    while (iss >> tmp_ceof) coef.emplace_back(stod(tmp_ceof));
    testPoly = Poly(coef);

    if (!getline(validfile, line)) break;
    istringstream iss2(line);
    while (iss2 >> tmp_ans) {
      if (tmp_ans == "#") {
        ans = {};
        break;
      } else {
        ans.emplace_back(stod(tmp_ans));
      }
    }
    cout << "Current Polynomial: " << testPoly << endl;
    if (validSinglePoly(testPoly, util, ans)) pass++;
    test_count++;
  }

  cout << "============================================"
       << " Total Test: " << test_count << " Passed Case: " << pass
       << " Pass Rate: " << pass / double(test_count) << endl;
  return;
}

void testBudan() {
  srand(time(NULL));
  validPolyFromFile("test/validation.test");

  return;
}

#endif