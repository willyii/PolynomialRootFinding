#ifndef VINCENT_H
#define VINCENT_H

#include "isolate.h"
#include <set>

using std::set;

class Vincent : public Isolate {
 public:
  vector<tuple<double, double>> isoRoot(Poly& p);
  vector<double> solve(Poly& p);
  Poly inverseTimesB(Poly& p, int b); // util

private:
};
#endif