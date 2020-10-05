#include "budan.h"

#include <unordered_map>

using namespace std;

// Poly after add h to x
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

int Budan::signChangeNum(Poly& p) {
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
