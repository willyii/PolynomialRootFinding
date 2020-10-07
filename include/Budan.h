#ifndef BUDAN_H
#define BUDAN_H

#include "isolate.h"

class Budan: public Isolate{
    public:
    int signChangeNum(Poly& tmp, double h);
    vector<double> solve(Poly& p);

    private:
    double bound(Poly& p);
    Poly addToP(Poly& p, double h);
    double rootInBound(Poly& p, double left, double right);
};

struct Boundry{
    double left;
    int lchange = 0; // sign change at left 
    double right;
    int rchange = 0; // sign change at right
};

#endif
