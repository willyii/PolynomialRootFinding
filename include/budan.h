#ifndef BUDAN_H
#define BUDAN_H

#include "isolate.h"

class Budan: public Isolate{
    public:
    vector<tuple<double, double>> isoRoot(Poly& p);
    vector<double> solve(Poly& p);
};


#endif
