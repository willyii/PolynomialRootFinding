//#include "bisect.h"

// vector<tuple<double, double>> bisectIsoroot(Poly &p) {
//  if (DEBUG_BISECT) cout << "DEBUG BISECT:    p: " << p << endl;
//  vector<tuple<double, double>> ret;
//  deque<tuple<double, double, Poly>> L = {std::make_tuple(0.0, 0.0, p)};
//  int N = p.deg();

//  while (L.size() > 0) {
//    double c = std::get<0>(L[0]);
//    double k = std::get<1>(L[0]);
//    Poly q = std::get<2>(L[0]);
//    if (DEBUG_BISECT) cout << "DEBUG BISECT:    c: " << c << endl;
//    if (DEBUG_BISECT) cout << "DEBUG BISECT:    k: " << k << endl;
//    if (DEBUG_BISECT) cout << "DEBUG BISECT:    q: " << q << endl;
//    L.pop_front();
//    if (q.valueAt(0) == 0) {
//      vector<double> tmpCoef = {1, 0};
//      Poly tt = Poly(tmpCoef);
//      q /= tt;
//      N -= 1;
//      ret.emplace_back(
//          std::make_tuple(c / pow(2, int(k)), (c + 0) / pow(2, int(k))));
//    }

//    Poly tmp = q;
//    tmp = reverseCoef(tmp);
//    tmp = addToP(tmp, 1);
//    if (DEBUG_BISECT) cout << "DDEBUG BISECT: tmp: " << tmp << endl;
//    sleep(1);
//    int v = signChangeNum(tmp, 0);

//    if (v == 1) {
//      ret.emplace_back(
//          std::make_tuple(c / pow(2, int(k)), (c + 1) / pow(2, int(k))));
//    }
//    if (v > 1) {
//      Poly tmpq = q;
//      tmpq = _change1(tmpq);
//      L.emplace_back(std::make_tuple(2 * c, k + 1, tmpq));

//      tmpq = addToP(tmpq, 1);
//      L.emplace_back(std::make_tuple(2 * c + 1, k + 1, tmpq));
//    }
//  }
//  return ret;
//}

// vector<double> bisectSolve(Poly &p) {
//  if (p.deg() == 1) return {-p[1] / p[0]};
//  vector<Poly> plist = squareFreeDecompo(p);
//  vector<tuple<double, double>> ranges, tmprange;
//  vector<double> roots;
//  double tmp;

//  if (DEBUG_BISECT) {
//    cout << "DEBUG BISECT SQUREFREE DECOMP: " << endl;
//    for (auto x : plist) cout << x << endl;
//  }

//  for (auto p : plist) {
//    if (p.deg() <= 0) continue;
//    tmprange = bisectIsoroot(p);
//    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
//  }

//  for (auto range : ranges) {
//    if (DEBUG_BISECT)
//      cout << "Finding root in range: " << get<0>(range) << " to "
//           << get<1>(range) << endl;
//    tmp = rootInBound(p, get<0>(range), get<1>(range));
//    if (tmp != NOTFOUND) roots.emplace_back(tmp);
//  }

//  return roots;
//}

// Poly reverseCoef(Poly &p) {
//  vector<double> coef = p.getCoef();
//  std::reverse(coef.begin(), coef.end());
//  Poly ret = Poly(coef);
//  return ret;
//}

// Poly _change1(Poly &p) {
//  vector<double> coef = p.getCoef();
//  double timer = 1;
//  int N = p.deg();
//  for (int i = 0; i <= N; i++) {
//    coef[i] *= timer;
//    timer *= 2;
//  }
//  Poly ret = Poly(coef);
//  return ret;
//}
