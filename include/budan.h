#ifndef BUDAN_H
#define BUDAN_H

#include "isolate.h"

class Budan: public Isolate{
    public:
    vector<double> solve(Poly& p);
};


#endif
