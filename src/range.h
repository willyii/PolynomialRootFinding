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

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/interval.hpp>
#include <iostream>

typedef boost::numeric::interval<double> interval;

struct Range {
  interval left_end;
  interval right_end;
};

/**
 * Print out Range
 */
template <int n>
static std::ostream &operator<<(std::ostream &out, const Range &u) {
  out << "left : " << u.left_end.lower()
      << " to right : " << u.right_end.upper();
  return out;
}

#endif // POLY_RANGE_H
