#ifndef SPARSE_H
#define SPARSE_H

#include <string>
#include <vector>

#include "poly.h"

struct Item {
  double c;  // coefficent
  int d;     // degree of x
};

std::vector<std::vector<double>> parseFromFile(char const* filename);

std::vector<double> parseSingleLine(std::string polyString);

Item parseItem(std::string item);

#endif
