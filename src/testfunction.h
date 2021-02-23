// ----------------------------------------------------------------------------
// FILENAME: testfunction.h
//
// DESCRIPTION:
//    This file defines and implement functions used to test the performance of
//    this project.
//
// FUNCTIONS:
//  random_int      :Get random integer from min to max
//  random_double   :Get random double from min to max
//  RandomPoly      :Get a polynomial with random coefficients
//  PolyToString    :Get string of this polynomial
//  RandomPolyToFile:Generate random polynomials and write to file kRANDOM_FILES
//  RunningTimeTest :Read polynomials from file and test the solving time of
//                   two methods
//
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#include "budan.h"
#include "parse.h"
#include "poly.h"
#include "range.h"
#include "util.h"
#include "vincent.h"

static const int digit = 1; // number of digit after point
static const int digit_control = std::pow(10, digit); // controler of digit
static const double kMAXROOT = 10;     // root from [-kMAXROOT, kMAXROOT]
static const double kCHANGE_PROB = .5; // prob to change root
static const bool FIXSEED = false;     // fix test set
static const int SEED = 2033;

/**
 * Get random integer in range min to max
 */
int rand_int(int min, int max) { return rand() % (max - min) + min; }

/**
 * Get random double in range min to max
 */
double rand_double(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  f = min + f * (max - min);
  f = std::ceil(f * digit_control) / digit_control;
  return f;
}

/**
 * struct RandomPolyRet - Return of RandomPoly
 */
struct RandomPolyRet {
  Poly<kMAXDEGREE> poly; // Polynomial
  vector<double> roots;  // Coressponding roots
};

/**
 * Generate a random polynomial.
 *
 */
RandomPolyRet RandomPoly() {
  RandomPolyRet ret;
  double root = 0.0;

  int num_roots = rand_int(2, kMAXDEGREE);
  while (root == 0.0)
    root = rand_double(-kMAXROOT, kMAXROOT);
  ret.poly[1] = interval(1, 1);
  ret.poly[0] = interval(-root, -root);
  ret.poly.set_degree(1);
  ret.roots.emplace_back(root);

  for (int i = 1; i < num_roots; i++) {
    if (rand_double(0, 1) < kCHANGE_PROB) {
      root = 0.0;
      while (root == 0.0)
        root = rand_double(-kMAXROOT, kMAXROOT);
    }

    Poly<kMAXDEGREE> backup(ret.poly);
    // Right shift ans by 1
    for (int j = ret.poly.get_degree(); j >= 0; j--)
      ret.poly[j + 1] = ret.poly[j];
    ret.poly[0] = 0;
    ret.poly.set_degree(ret.poly.get_degree() + 1);
    ret.poly -= backup * root;

    ret.roots.emplace_back(root);
  }

  return ret;
}

/**
 * Return corresponding string of polynomial
 */
template <int n> std::string PolyToString(const Poly<n> &poly) {
  std::string ret = "";
  for (int i = 0; i <= n; i++) {
    if (poly.containZero(i))
      continue;
    if (poly[i] >= 0)
      ret += '+';
    ret += std::to_string(boost::numeric::median(poly[i])) + "*x^" +
           std::to_string(i);
  }
  return ret;
}

/**
 * Generate some polynomials and save to file
 *
 * @param num_polys :Number of polynomials to generate
 */
void RandomPolyToFile(int num_polys) {

  srand(time(NULL));

  if (FIXSEED)
    srand(SEED);
  std::ofstream write_file_poly, write_file_solution;
  write_file_poly.open(kRANDOM_FILE);
  write_file_solution.open(kRANDOM_FILE_SOL);

  for (int i = 0; i < num_polys; i++) {
    auto single_poly = RandomPoly();
    write_file_poly << PolyToString(single_poly.poly) << "\n";

    std::string root_string("");
    for (auto root : single_poly.roots)
      root_string += std::to_string(root) + ",";

    write_file_solution << root_string.substr(0, root_string.length() - 1)
                        << "\n";
  }

  write_file_poly.close();
  write_file_solution.close();
  return;
}

/**
 * Get polynomials from file, solve use budan't theorem and continued fractions.
 * print running time of these two methods
 *
 * @param file_path :Path to file
 */
void RunningTimeTest(const char *file_path) {
  vector<double *> coefs;
  vector<int> num_coefs;
  vector<Range *> budan_roots;
  vector<int> budan_root_num;
  vector<Range *> vincent_roots;
  vector<int> vincent_root_num;

  Range *poly_roots(nullptr);
  ParseFromFile(file_path, coefs, num_coefs);

  // check num of polynoals
  assert(coefs.size() == num_coefs.size());

  /// Budan's theorem
  auto budan_start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < coefs.size(); i++) {
    poly_roots = new Range[kMAXDEGREE];
    budan_root_num.push_back(
        BudanRootIsolate(coefs[i], num_coefs[i], poly_roots));
    budan_roots.push_back(poly_roots);
    delete[] poly_roots;
  }
  auto budan_end = std::chrono::high_resolution_clock::now();

  auto budan_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      budan_end - budan_start);

  /// continued fration
  auto vincent_start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < coefs.size(); i++) {
    poly_roots = new Range[kMAXDEGREE];
    vincent_root_num.push_back(
        VincentRootIsolate(coefs[i], num_coefs[i], poly_roots));
    vincent_roots.push_back(poly_roots);
    delete[] poly_roots;
  }
  auto vincent_end = std::chrono::high_resolution_clock::now();

  auto vincent_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      vincent_end - vincent_start);

  // Print budan time
  std::cout << "Buand Method takes " << budan_duration.count() << " ms for "
            << coefs.size() << " polynomials" << std::endl;
  std::cout << "Average time : " << budan_duration.count() / coefs.size()
            << "ms" << std::endl;

  // Print continued fraction time
  std::cout << "Vincent Method takes " << vincent_duration.count() << " ms "
            << coefs.size() << " polynomials" << std::endl;
  std::cout << "Average time : " << vincent_duration.count() / coefs.size()
            << "ms" << std::endl;
  return;
}

/**
 * @brief Read roots from file
 */
vector<vector<double>> ParseRootFile(const char *root_file_path) {
  string lineString;
  std::ifstream fileString(root_file_path);
  vector<vector<double>> ret;

  while (getline(fileString, lineString)) { // every line
    vector<double> current_root;

    std::stringstream ss(lineString);
    string root_string;
    while (getline(ss, root_string, ','))
      current_root.emplace_back(stod(root_string));

    ret.emplace_back(current_root);
  }

  return ret;
}

/**
 * Get polynomials from file, solve use budan't theorem and continued
 * fractions. print correct rate of these two methods
 *
 * @param file_path :Path to file
 */
void RunningCorrectTest(const char *test_file_path,
                        const char *valid_file_path = "") {
  vector<double *> coefs;
  vector<int> num_coefs;
  vector<vector<double>> actual_roots;
  int root_num;

  // Read polynomials from file
  Range *poly_roots(nullptr);
  ParseFromFile(test_file_path, coefs, num_coefs);

  // Read roots from file
  if (std::string(valid_file_path) != "")
    actual_roots = ParseRootFile(valid_file_path);

  // check num of polynoals
  assert(coefs.size() == num_coefs.size());

  for (size_t i = 0; i < coefs.size(); i++) {
    std::cout << "=========== " << i << " =============" << std::endl;

    std::cout << " Polynomial: " << Poly<kMAXDEGREE>(coefs[i], num_coefs[i])
              << std::endl;
    if (std::string(valid_file_path) != "") {
      printf("actual roots: \n");
      for (double root : actual_roots[i])
        printf("%f , ", root);
      printf("\n \n");
    }

    /// ----------------------
    /// Budan theorem
    /// ----------------------
    poly_roots = new Range[kMAXDEGREE];
    root_num = BudanRootIsolate(coefs[i], num_coefs[i], poly_roots);
    std::cout << "root num " << root_num << std::endl;
    // assert(root_num == actual_roots[i].size());
    printf("Budan Theorem resuts: \n");
    for (int j = 0; j < root_num; j++)
      printf("From %f to %f \n", boost::numeric::median(poly_roots[j].left_end),
             boost::numeric::median(poly_roots[j].right_end));
    delete[] poly_roots;
    poly_roots = nullptr;
    printf("\n");

    /// ----------------------
    /// Continued fractions
    /// ----------------------
    poly_roots = new Range[kMAXDEGREE];
    root_num = VincentRootIsolate(coefs[i], num_coefs[i], poly_roots);
    printf("Continued Fraction resuts: \n");
    for (int j = 0; j < root_num; j++)
      printf("From %f to %f \n", boost::numeric::median(poly_roots[j].left_end),
             boost::numeric::median(poly_roots[j].right_end));
    delete[] poly_roots;
    printf("\n \n \n");
  }
  return;
}
