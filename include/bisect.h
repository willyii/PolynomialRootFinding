#ifndef BISECT_H
#define BISECT_H

#include "util.h"

vector<tuple<double, double>> bisectIsoroot(Poly &p);

vector<double> bisectSolve(Poly &p);

Poly reverseCoef(Poly &p);

Poly _change1(Poly &p);

#endif
