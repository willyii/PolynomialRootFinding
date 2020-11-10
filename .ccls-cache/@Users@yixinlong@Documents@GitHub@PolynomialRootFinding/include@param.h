#ifndef PARAM_H
#define PARAM_H
#include <unistd.h>

#include <cstdlib>
#include <limits>

// #define EPSILON std::numeric_limits<double>::epsilon()

#define NOTFOUND std::numeric_limits<double>::min()

#define EPSILON 1e-9
#define TESTERROR 1e-5
#define MAXITER 100000
#define MINRANGE 0.000001

#define DEBUG_GCD false
#define DEBUG_SQD false
#define DEBUG_BUDAN false
#define DEBUG_VINCENT true
#define DEBUG_BISECT false

#endif