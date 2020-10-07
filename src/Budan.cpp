#include "budan.h"

#include <deque>
#include <unordered_map>

using namespace std;

Poly Budan::addToP(Poly& p, double h) {
  int N = p.size();
  Coef tmp = Coef(N);
  unordered_map<int, Coef> memo;  // map the idx to corresponding coef

  // zero case
  memo[0] = Coef(N);
  memo[0][N - 1] = 1;

  // one case
  if (N >= 1) {
    memo[1] = Coef(N);
    memo[1][N - 1] = h;
    memo[1][N - 2] = 1;
  }

  // more case
  for (int idx = 2; idx < N; idx++) {
    int mid = idx / 2;
    memo[idx] = memo[mid] * memo[idx - mid];
  }

  // combine
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      tmp[j] += p[N - i - 1] * memo[i][j];
    }
  }

  return Poly(tmp);
}

int Budan::signChangeNum(Poly& tmp, double h) {
  Poly p = addToP(tmp, h);
  int ret = 0;
  bool sign = p[0] > 0;
  for (int i = 0; i < p.size(); i++) {
    if (abs(p[i] - 0) <= EPSILON) continue;
    if (p[i] > 0 == sign) continue;
    sign = !sign;
    ret++;
  }
  return ret;
}

// TODO: solve function
vector<double> Budan::solve(Poly& p) {
  vector<double> roots;
  deque<Boundry> b;
  double tmp = bound(p), mid;
  int midchange;

  b.push_back(
      Boundry{-tmp, signChangeNum(p, -tmp), tmp, signChangeNum(p, tmp)});
  Boundry tmpb;

  // TODO: Apply bisection to find the range of the root
  while (b.size() > 0) {
    tmpb = b[0];
    b.pop_front();
    mid = (tmpb.left + tmpb.right) / 2;
    midchange = signChangeNum(p, mid);

    // left side
    if (tmpb.lchange - midchange == 1) {
      roots.push_back(rootInBound(p, tmpb.left, mid));
      // roots.push_back((tmpb.left + mid) / 2);  // TODO solve exactly
    } else if (tmpb.lchange - midchange > 1) {
      b.push_back(Boundry{tmpb.left, tmpb.lchange, mid, midchange});
    }

    // right side
    if (midchange - tmpb.rchange == 1) {
      roots.push_back(rootInBound(p, mid, tmpb.right));
      // roots.push_back((tmpb.right + mid) / 2);  // TODO solve exactly
    } else if (midchange - tmpb.rchange > 1) {
      b.push_back(Boundry{mid, midchange, tmpb.right, tmpb.rchange});
    }
  }

  return roots;
}

// Get the boundry of roots, applied Cauchy's bound
double Budan::bound(Poly& p) {
  int N = p.size();
  double tmp = __DBL_MIN__, lc = p.getCoef().lc();
  for (int i = 0; i < N - 1; i++) {
    tmp = max(tmp, fabs(p[i] / lc));
  }
  return 1 + tmp;
}

// Finding root in boundry b. root is isolated in that boundry
// Newtown method temporary
double Budan::rootInBound(Poly& p, double left, double right){
  double x0 = (left + right)/2;
  while(p.valueAt(x0) != 0 ){
    x0 = x0- (p.valueAt(x0)/p.gradientAt(x0));
  }
  return x0;
}
