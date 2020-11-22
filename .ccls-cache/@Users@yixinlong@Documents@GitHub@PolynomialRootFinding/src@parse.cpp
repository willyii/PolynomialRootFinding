#include "parse.h"

#include <array>
#include <fstream>
#include <iostream>

/*
 * Read from file, format like x^4-1.2*x^3+.51*x^2-.88e-1*x+.48e-2 
 * which cotains many polys
 * this should return vector of coef,
 * this coef should not lead with 0, make no assumption of max degree
 */
std::vector<std::vector<double>> parseFromFile(char const* filename) {
  std::vector<std::vector<double>> coef_list;

  std::string polyString;
  std::ifstream fileString(filename);
  while (getline(fileString, polyString)) {
    /* strip space */
    polyString.erase(remove(polyString.begin(), polyString.end(), ' '),
                     polyString.end());
    std::cout << polyString << std::endl;

    // std::vector<double> coef = parseSingleLine(polyString);
    coef_list.emplace_back(parseSingleLine(polyString));
  }

  return coef_list;
}

/*
 * This function take single polynomial string as input
 * return the coef of this polynomail, constant last
 * this coefficent should not assume the max degree
 */
std::vector<double> parseSingleLine(std::string polyString) {
  std::string item = "";
  std::vector<Item> itemList;

  for (size_t i = 0; i < polyString.length(); i++) {
    if ((polyString[i] != '+' && polyString[i] != '-') ||
        polyString[i - 1] == 'e')
      item += polyString[i];
    else {
      itemList.push_back(parseItem(item));
      item = polyString[i];
    }
  }
  if (item.length() > 0) {
    itemList.push_back(parseItem(item));
  }

  // maxd
  int maxD = 0;
  for (auto i : itemList) maxD = std::max(maxD, i.d);
  // std::cout << "maxD: " << maxD << std::endl;
  vector<double> coef(maxD + 1, 0.0);
  for (auto i : itemList) coef[i.d] += i.c;

  std::reverse(coef.begin(), coef.end());
  // for (size_t i = 0; i < coef.size(); i++) {
  //  std::cout << "degree: " << i << ", coeff: " << coef[i] << std::endl;
  //}

  return coef;
}

/*
 * This function used for parse single item like in format c*x^d
 * And this should return Item struct with c = c and d = d
 */
Item parseItem(std::string item) {
  // std::cout << "item = " << item << std::endl;
  Item ans;  // save the result
  double c;  // coefficient of this item
  int d;     // degree of this tiem

  std::size_t pos = item.find('x');
  if (pos == std::string::npos) {  // Constant
    c = std::stod(item);
    d = 0;
  } else if (pos == (item.length() - 1)) {  // x
    c = std::stod(item.substr(0, item.length() - 2));
    d = 1;
  } else if (pos == 0) {  // x^d
    c = 1.0;
    d = std::stoi(item.substr(item.length() - 1));
  } else {  // c*x^d
    c = std::stod(item.substr(0, item.length() - 4));
    d = std::stoi(item.substr(item.length() - 1));
  }

  // std::cout << "d = " << d << ", c = " << c << std::endl;
  ans.d = d;
  ans.c = c;
  return ans;
}
