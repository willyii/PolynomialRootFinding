#ifndef BUDAN_H
#define BUDAN_H

#include "isolate.h"

class Budan: public Isolate{
    public:
    Poly addToP(Poly& p, double h);
    int signChangeNum(Poly& p);

};

#endif
