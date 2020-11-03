#ifndef BISECT_H
#define BISECT_H

#include "util.h"

vector<tuple<double, double>> bisectIsoroot(Poly &p);

vector<double> bisectSolve(Poly &p);

#endif
