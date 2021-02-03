// ----------------------------------------------------------------------------
//
// FILENAME: parse.h
//
// DESCRIPTION:
//    This file defines and implement the functions used to parse polynomials in
//    file
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------
#ifndef POLY_PARSE_H
#define POLY_PARSE_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "poly.h"

using std::string;
using std::vector;

/**
 * @brief Parse item string, eg: +1.23*x^5
 *
 * @param item_string :String of item, might be modified
 * @param coef :Pointer to coefficents
 * @return degree of this item
 */
int ParseItem(string &item_string, double *coef) {
  /// postion of x
  std::size_t pos = item_string.find('x');

  if (pos == string::npos) { // constant
    coef[0] += std::stod(item_string);
    return 0;
  }

  if (pos == item_string.length() - 1) // end with x
    item_string += "^1";

  if (pos == 0) /// x..
    item_string = "+1*" + item_string;
  else if (pos == 1) {
    if (item_string[0] == '+') /// +x ..
      item_string = "+1*" + item_string;
    else /// -x...
      item_string = "-1*" + item_string;
  }

  pos = item_string.find('x');
  int d = std::stoi(item_string.substr(pos + 2));
  coef[d] += std::stod(item_string.substr(0, pos));

  return d;
}

/**
 * @brief :Parse a line like 1.23*x^6 + 5*x +2
 *
 * @param polyString :String to parse
 * @param coef :Pointer to coefficents, might be modified
 * @return :Largest degree of this polynomial
 */
int ParseSingleLine(const string &polyString, double *coef) {
  string item_string("");
  int ret(0);

  for (int i = 0; i <= kMAXDEGREE; i++)
    coef[i] = 0.0;

  for (size_t i = 0; i < polyString.length(); i++) {
    if ((polyString[i] != '+' && polyString[i] != '-') || i == 0 ||
        polyString[i - 1] == 'e')
      item_string += polyString[i];
    else {
      // Process item_string
      ret = std::max(ParseItem(item_string, coef), ret);
      item_string = polyString[i];
    }
  }

  if (item_string.length() > 0)
    ret = std::max(ParseItem(item_string, coef), ret);

  return ret;
}

/**
 * @brief :Parse polynomials in file
 *
 * @param filename :Path to file
 * @param ceofs :vector of pointers that point to coefficents, might be modified
 * @param num_coefs :vector of int, ref to the number of coef coresponding to
 * coefs
 */
void ParseFromFile(const char *filename, vector<double *> &coefs,
                   vector<int> &num_coefs) {
  string polyString;
  double *coef(nullptr);
  int line_count = 1;

  // std::ifstream fileString("src/test.poly");
  std::ifstream fileString(filename);

  while (getline(fileString, polyString)) { /// every line

    /* strip space */
    polyString.erase(std::remove(polyString.begin(), polyString.end(), ' '),
                     polyString.end());

    /* Print out current polynomial */
    std::cout << line_count++ << ": " << polyString << std::endl;

    coef = new double[kMAXDEGREE + 1];

    num_coefs.push_back(ParseSingleLine(polyString, coef) + 1);

    coefs.push_back(coef);
  }

  return;
}

#endif // POLY_PARSE_H
