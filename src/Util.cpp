#include "Util.h"

#include <unordered_map>
#include <vector>

#include "Param.h"

using namespace std;
using namespace util;

// Check if Zero vector or not 
bool util::isZeroVec(vector<double> &c) {
  for (auto num : c) {
    if (abs(num) > Param::EPSILON) return false;
  }
  return true;
}

// Leading coeff of Polynomial
double util::lc(vector<double> &c) {
  for (auto num : c) {
    if (num != 0.0) {
      return num;
    }
  }
  return 0.0;
}

// Degree of Polynomial
int util::deg(vector<double> &c) {
  int ret = 0;
  for (int i = 0; i < c.size(); i++) {
    if (c[i] != 0.0) {
      return c.size() - i - 1;
    }
  }
  return 0;
}

// Caculate the coefficient after add h to x
vector<double> util::coeffAfter(vector<double> &coef, double h) {
  int N = coef.size();
  vector<double> tmp(N, 0);
  unordered_map<int, vector<double>> memo;  // map the idx to corresponding coef

  // zero case
  memo[0] = vector<double>(N, 0);
  memo[0][N - 1] = 1;

  // one case
  if (N >= 1) {
    memo[1] = vector<double>(N, 0);
    memo[1][N - 1] = h;
    memo[1][N - 2] = 1;
  }

  // more case
  for (int idx = 2; idx < N; idx++) {
    int mid = idx / 2;
    memo[idx] = polyTimes(memo[mid], memo[idx - mid]);
  }

  // combine
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      tmp[j] += coef[N - i - 1] * memo[i][j];
    }
  }

  return tmp;
}

// Calculate the coef after two poly times
vector<double> util::polyTimes(vector<double> &c1, vector<double> &c2) {
  vector<double> ret(c1.size(), 0);
  int idx1, idx2;
  for (int i = 0; i < c1.size(); i++) {
    idx1 = c1.size() - i - 1;
    for (int j = 0; j < c1.size(); j++) {
      idx2 = c1.size() - j - 1;
      if (idx1 + idx2 >= c1.size()) continue;
      ret[c1.size() - 1 - idx1 - idx2] += c1[i] * c2[j];
    }
  }
  return ret;
}

// Calculate the sum of two poly coef
vector<double> util::polyAdd(vector<double> &c1, vector<double> &c2) {
  vector<double> ret(c1.size(), 0);
  for (int i = 0; i < c1.size(); i++) {
    ret[i] = c1[i] + c2[i];
  }
  return ret;
}

// Calculate the sub of two poly coef
vector<double> util::polySub(vector<double> &c1, vector<double> &c2) {
  vector<double> ret(c1.size(), 0);
  for (int i = 0; i < c1.size(); i++) {
    ret[i] = c1[i] - c2[i];
  }
  return ret;
}

// Calculate the div of two poly coef
vector<double> util::polyDiv(vector<double> &c1, vector<double> &c2) {
    vector<double> q(c1.size(), 0.0), r = c1,t, tmp;
    while(!isZeroVec(r) and deg(r) >= deg(c2)){
        t = leadDiv(r, c2);
        q = polyAdd(r,t);
        tmp = polyTimes(t,c2);
        r = polySub(r, tmp);
    }
    return q;
}

// Divide the leading term
vector<double> util::leadDiv(vector<double> &c1, vector<double> &c2) {
    double div = lc(c1)/ lc(c2);
    int degree = deg(c1) - deg(c2);
    vector<double> t (c1.size(), 0.0);
    t[c1.size() - degree - 1] = div;
    return t;
}