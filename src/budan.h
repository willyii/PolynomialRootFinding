#ifndef BUDAN_H
#define BUDAN_H

#include "boundry.h"
#include "util.h"

using std::tuple;

// Isolate the range of every root
vector<tuple<double, double>> budanIsoroot(Poly& p);

// Slover using Budan's theorem
vector<double> budanSolve(Poly& p);

#endif
