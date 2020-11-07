#ifndef VINCENT_H
#define VINCENT_H

#include "util.h"

using std::tuple;

// Isolat the range of every root
vector<tuple<double, double>> vincentIsoroot(Poly& p, bool negative);

// Solver using Vincent
vector<double> vincentSolve(Poly& p);

// Inverse the coefficent and generate new Poly
Poly reverse(Poly& p);

// Refine the range, make the range smaller enough to avoid edge case
void refineRange(Poly& p, tuple<double, double>& range);

// Interval
struct Interval {
  double a, b, c, d;  // Set range
  Poly p;
  int s;  // Sign Change
};

#endif
