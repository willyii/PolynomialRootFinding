// ----------------------------------------------------------------------------
//
// FILENAME: range.h
//
// DESCRIPTION:
//    This file define the struct to store the range of root
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_RANGE_H
#define POLY_RANGE_H

#include <iostream>

struct Range {
  double left_end;
  double right_end;
};

/**
 * Print out Range
 */
template <int n>
static std::ostream &operator<<(std::ostream &out, const Range &u) {
  out << "left : " << u.left_end << " to right : " << u.right_end;
  return out;
}

#endif // POLY_RANGE_H
