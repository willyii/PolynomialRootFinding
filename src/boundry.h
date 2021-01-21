#ifndef POLY_BOUNDRY_H
#define POLY_BOUNDRY_H

struct Boundry {
  double left;
  int lchange;  // sign change at left
  double right;
  int rchange;  // sign change at right
};

#endif  // POLY_BOUNDRY_H
