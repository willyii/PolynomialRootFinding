#include "budan.h"

#include <deque>
#include <unordered_map>

using namespace std;

Poly Budan::addToP(Poly& p, double h) {
  int N = p.size();
  vector<double> tmp(N, 0.0);
  unordered_map<int, Poly> memo;  // map the idx to corresponding coef
  // zero case
  tmp[N - 1] = 1;
  memo[0] = Poly(tmp);
  fill(tmp.begin(), tmp.end(), 0);
  // one case
  if (N >= 1) {
    tmp[N - 1] = h;
    tmp[N - 2] = 1;
    memo[1] = Poly(tmp);
    fill(tmp.begin(), tmp.end(), 0);
  }
  // more case
  for (int idx = 2; idx < N; idx++) {
    int mid = idx / 2;
    memo[idx] = memo[mid] * memo[idx - mid];
  }
  // combine
  Poly ans = Poly(tmp);
  for (int i = 0; i < N; i++) {
    memo[i] = memo[i] * p[N - i - 1];
    ans += memo[i];
  }
  if (DEBUG_BUDAN)
    cout << "DEBUG BUDAN: after add " << h << " : " << Poly(tmp) << endl;
  return ans;
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

vector<tuple<double, double>> Budan::isoRoot(Poly& p) {
  vector<tuple<double, double>> ret;
  deque<Boundry> b;
  double tmp = bound(p), mid, root;
  int midchange;
  Boundry tmpb;

  b.push_back(
      Boundry{-tmp, signChangeNum(p, -tmp), tmp, signChangeNum(p, tmp)});

  while (b.size() > 0) {
    tmpb = b[0];
    b.pop_front();
    if (DEBUG_BUDAN)
      cout << "Search Boundary: " << tmpb.left << "\t | \t  " << tmpb.right
           << endl;
    if (DEBUG_BUDAN)
      cout << tmpb.lchange << "\t | \t  " << tmpb.rchange << endl;
    mid = (tmpb.left + tmpb.right) / 2;
    midchange = signChangeNum(p, mid);
    if (DEBUG_BUDAN) cout << "Mid Change:  " << midchange << endl;

    // left side
    if (mid - tmpb.left < MINRANGE && tmpb.lchange - midchange > 0) {
      root = rootInBound(p, tmpb.left, mid);
      if (root != NOTFOUND)
        ret.push_back(tuple<double, double>(tmpb.left, mid));
    } else if (tmpb.lchange - midchange > 0) {
      b.push_back(Boundry{tmpb.left, tmpb.lchange, mid, midchange});
    }

    // right side
    if (tmpb.right - mid < MINRANGE && midchange - tmpb.rchange > 0) {
      root = rootInBound(p, mid, tmpb.right);
      if (root != NOTFOUND)
        ret.push_back(tuple<double, double>(mid, tmpb.right));
    } else if (midchange - tmpb.rchange > 0) {
      b.push_back(Boundry{mid, midchange, tmpb.right, tmpb.rchange});
    }
  }

  return ret;
}

// Get the boundry of roots, applied Cauchy's bound
double Budan::bound(Poly& p) {
  int N = p.size();
  double tmp = __DBL_MIN__, lc = p.lc();
  for (int i = 1; i < N; i++) {
    tmp = fmax(tmp, fabs(p[i] / lc));
  }
  return 1 + tmp;
}

// Finding root in boundry b. root is isolated in that boundry
// Newtown method temporary
double Budan::rootInBound(Poly& p, double left, double right) {
  double x0 = (left + right) / 2;
  int idx = 0;
  while (p.valueAt(x0) != 0 && idx < MAXITER && x0 >= left && x0 <= right) {
    x0 = x0 - (p.valueAt(x0) / p.gradientAt(x0));
    idx += 1;
  }
  if (idx == MAXITER || x0 < left || x0 > right) {
    return NOTFOUND;
  }
  return x0;
}

vector<double> Budan::solve(Poly& p) {
  vector<Poly> plist = squareFreeDecompo(p);
  vector<tuple<double, double>> ranges, tmprange;
  vector<double> roots;
  double tmp;

  if (DEBUG_BUDAN) {
    cout << "DEBUG BUDAN SQUREFREE DECOMP: " << endl;
    for (auto x : plist) cout << x << endl;
  }

  for (auto p : plist) {
    if (p.deg() <= 0) continue;
    tmprange = isoRoot(p);
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_BUDAN)
      cout << "Finding root in range: " << get<0>(range) << " to "
           << get<1>(range) << endl;
    tmp = rootInBound(p, get<0>(range), get<1>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }

  if (roots.size() == 0)
    cout << "There is no real roots of this polynomial" << endl;
  return roots;
}
