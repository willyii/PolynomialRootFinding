#include "Budan.h"

#include <unordered_map>
#include <vector>
#include <iostream>

#include "Param.h"

using namespace std;

/**
 * Caculate the number of sign change when this equation added h
 *
 * @param poly Polynomila function that need to be caculated
 * @param h after add this number, formula will be changed
 * @return the number of sign change
 */
int Budan::signChangeNums(Poly &poly, double h) {
  int ret = 0;
  vector<double> after = coeffAfter(poly, h);

  bool sign = after[0] > 0;
  for (auto num : after) {
    if (num == 0) continue;
    if (num > 0 == sign) continue;
    sign = !sign;
    ret++;
  }
  return ret;
}

/**
 * Caculate the coefficient after add h to x
 *
 * @param poly Polynomial function
 * @param h the change of x
 * @return coef after change
 */
vector<double> Budan::coeffAfter(Poly &poly, double h) {
  int N = poly.getN();
  vector<double> tmp(N, 0), coef = poly.getCoef();
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

/**
 * Calculate the coef after two poly times
 * We can make sure the maximum will no large than _N
 *
 * @param c1 coef of poly1
 * @param c2 coef of poly2
 * @return the cofe after computation
 */
vector<double> Budan::polyTimes(vector<double> &c1, vector<double> &c2) {
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

/**
 * NewtonRaphson method to find root.
 *
 * This is already garentee that only one root in polynomoal.
 * TODO: Add intervals to this function and make sure it will find the root.
 *
 * @param poly polynomial function
 * @param x initial guess of the root.
 * @return the root of this polynomial.
 */
double Budan::NewtonRaphson(Poly &poly, double x) {
  int idx = 0;
  while ((abs(poly.getValue(x) - 0) > Param::EPSILON) &&
         (idx < Param::MAX_ITER)) {
    x -= poly.getValue(x) / poly.getGradient(x);
    idx += 1;
  }
  if (abs(poly.getValue(x) - 0) <= Param::EPSILON) return x;
  return numeric_limits<double>::max();
}

/**
 * Finding root
 * 
 * @param poly
 * @return roots
 */
vector<double> Budan::rootFinding(Poly& poly){
  return vector<double> {0.0};
}


/**
 * Finding GCD of two Poly's coefficient by Euclidean division
 *
 * @param c1, c2, coeff of two polys
 * @return gcd of two polys
 */
GCD Budan::gcd(vector<double> &a, vector<double> &b){
  int N = a.size();
  vector<double> q(N, 0.0), r = a; 
  int d = Budan::deg(b);
  double c = Budan::lc(b);

  while(Budan::deg(r) >= d){
    vector<double> s(N, 0.0);
    s[N - deg(r) + d - 1] = Budan::lc(r)/c;
    vector<double> sb = Budan::polyTimes(s, b);
    for(int i=0;i<N;i++){
      q[i] = q[i] + s[i];
      r[i] = r[i] - sb[i];
    }
  }
  GCD tmp;
  tmp.q = q;
  tmp.r = r;
  return tmp; 
} 


double Budan::lc(vector<double> &c){
  for(auto num:c){
    if(num!= 0.0){
      return num;
    }
  }
  return 0.0;
}

int Budan::deg(vector<double> &c){
  int ret = 0;
  for(int i=0; i<c.size();i++){
    if(c[i] != 0.0){
      return c.size() - i -1;
    }
  }
  return 0;
}





