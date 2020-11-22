#include "parse.h"

#include <array>
#include <fstream>
#include <iostream>

std::vector<double*> parseFromFile(char const* filename) {
  std::vector<double*> coef_list;

  std::string polyString;
  std::ifstream fileString(filename);
  while (getline(fileString, polyString)) {
    /* strip space */
    polyString.erase(remove(polyString.begin(), polyString.end(), ' '),
                     polyString.end());
    std::cout << polyString << std::endl;

    // double coef[] = {0, 0, 0, 0, 0, 0, 0};  // Assume degree is 6
    double* coef = new double[7];
    for (size_t i = 0; i < 7; i++) coef[i] = 0.0;
    parseSingleLine(coef, polyString);
    coef_list.emplace_back(std::move(coef));
  }

  return coef_list;
}

void parseSingleLine(double* coef, std::string polyString) {
  // for (size_t i = 0; i < 7; i++) coef[i] = 0.0;
  std::string item = "";

  for (size_t i = 0; i < polyString.length(); i++) {
    if ((polyString[i] != '+' && polyString[i] != '-') ||
        polyString[i - 1] == 'e')
      item += polyString[i];
    else {
      parseItem(coef, item);
      item = polyString[i];
    }
  }
  if (item.length() > 0) parseItem(coef, item);

  for (size_t i = 0; i < 7; i++) {
    std::cout << "degree: " << i << ", coeff: " << coef[i] << std::endl;
  }

  return;
}

void parseItem(double* coef, std::string item) {
  std::cout << "item = " << item << std::endl;
  double c;  // coefficient of this item
  int d;     // degree of this tiem

  std::size_t pos = item.find('x');
  if (pos == std::string::npos) {  // Constant
    c = std::stod(item);
    d = 6;
  } else if (pos == (item.length() - 1)) {  // x
    c = std::stod(item.substr(0, item.length() - 1));
    d = 5;
  } else if (pos == 0) {  // x^d
    c = 1.0;
    d = 7 - std::stoi(item.substr(item.length() - 1)) - 1;
  } else {  // c*x^d
    d = 7 - std::stoi(item.substr(item.length() - 1)) - 1;
    c = std::stod(item.substr(0, item.length() - 4));
  }

  std::cout << "d = " << d << ", c = " << c << std::endl;
  for (size_t i = 0; i < d; i++) coef++;
  *coef = c;
}
