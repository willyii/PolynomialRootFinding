#ifndef BOUNDRY_H
#define BOUNDRY_H
struct Boundry {
  double left;
  int lchange;  // sign change at left
  double right;
  int rchange;  // sign change at right
};
#endif