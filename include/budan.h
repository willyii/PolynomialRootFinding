#ifndef BUDAN_H
#define BUDAN_H

#include "isolate.h"
#include <tuple>

using std::tuple;

class Budan: public Isolate{
    public:
    int signChangeNum(Poly& tmp, double h);
    vector<tuple<double, double>> isoRoot(Poly& p);
    vector<double> solve(Poly& p);

    private:
    double bound(Poly& p);
    Poly addToP(Poly& p, double h);
    double rootInBound(Poly& p, double left, double right);
};

struct Boundry{
    double left;
    int lchange; // sign change at left 
    double right;
    int rchange; // sign change at right
};

#endif
