#ifndef UTIL_H
#define UTIL_H
#include <vector>

using std::vector;

namespace util{
    bool isZeroVec(vector<double> &c);
    double lc(vector<double> &c);
    int deg(vector<double> &c);
    vector<double> coeffAfter(vector<double> &coef, double h);
    vector<double> polyTimes(vector<double> &c1, vector<double> &c2);
    vector<double> polyAdd(vector<double> &c1, vector<double> &c2);
    vector<double> polySub(vector<double> &c1, vector<double> &c2);
    vector<double> polyDiv(vector<double> &c1, vector<double> &c2);

    // util for util
    vector<double> leadDiv(vector<double> &c1, vector<double> &c2);
};

#endif