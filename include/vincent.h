#ifndef VINCENT_H
#define VINCENT_H

#include "isolate.h"

class Vincent : public Isolate {
 public:
  vector<tuple<double, double>> isoRoot(Poly& p);
  vector<double> solve(Poly& p);
};
#endif