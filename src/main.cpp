#include "budan.h"
#include "parse.h"
#include "poly.h"
#include "randompoly.h"
#include "range.h"
#include "util.h"
#include "vincent.h"

#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <vector>

/**
 * Get polynomials from file, solve use budan't theorem and continued fractions.
 * print running time of these two methods
 *
 * @param file_path :Path to file
 */
void RunningTime(const char *file_path) {
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
  }
  auto budan_end = std::chrono::high_resolution_clock::now();

  auto budan_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      budan_end - budan_start);

  // Print budan time
  std::cout << "Buand Method takes " << budan_duration.count() << " ms for "
            << coefs.size() << " polynomials" << std::endl;
  std::cout << "Average time : " << budan_duration.count() / coefs.size()
            << "ms" << std::endl;

  /// continued fration
  auto vincent_start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < coefs.size(); i++) {
    poly_roots = new Range[kMAXDEGREE];
    vincent_root_num.push_back(
        VincentRootIsolate(coefs[i], num_coefs[i], poly_roots));
    vincent_roots.push_back(poly_roots);
  }
  auto vincent_end = std::chrono::high_resolution_clock::now();

  auto vincent_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      vincent_end - vincent_start);

  // Print continued fraction time
  std::cout << "Vincent Method takes " << vincent_duration.count() << " ms "
            << coefs.size() << " polynomials" << std::endl;
  std::cout << "Average time : " << vincent_duration.count() / coefs.size()
            << "ms" << std::endl;
  return;
}

int main(int argc, char *argv[]) {
  if (argc == 1) { // get random polys
    for (int i = 0; i < 100; i++) {
      RandomPolyToFile();
      RunningTime(kRANDOM_FILE);
    }
  } else if (argc == 2) {
    RunningTime(argv[1]);
  }

  return 0;
}
