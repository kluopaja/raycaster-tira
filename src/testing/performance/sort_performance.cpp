#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <string>
#include <random>

#include "../../utils.h"
#include "../../vector.h"
#include "performance_utils.h"
#include "../../test_utils.h"

namespace {
struct SortTestResult {
  double qsort_time = 0.0;
  double std_sort_time = 0.0;
  double quick_sort_time = 0.0;
  double mm_quick_sort_time = 0.0;
};
SortTestResult operator+(const SortTestResult& a, const SortTestResult& b) {
  return {a.qsort_time + b.qsort_time, a.std_sort_time + b.std_sort_time,
          a.quick_sort_time + b.quick_sort_time,
          a.mm_quick_sort_time + b.mm_quick_sort_time};
}
SortTestResult operator/(const SortTestResult& a, double b) {
  return {a.qsort_time / b, a.std_sort_time / b, a.quick_sort_time / b,
          a.mm_quick_sort_time / b};
}
std::ostream& operator<<(std::ostream& out, const SortTestResult& result) {
  out << std::left;
  out << std::setw(25);
  out << "qsort: " << result.qsort_time/1e6 << "ms" << std::endl;
  out << std::setw(25);
  out << "std::sort: " << result.std_sort_time/1e6 << "ms" << std::endl;;
  out << std::setw(25);
  out << "quickSort: " << result.quick_sort_time/1e6 << "ms" << std::endl;
  out << std::setw(25);
  out << "mmQuickSort: " << result.mm_quick_sort_time/1e6 << "ms" << std::endl;
  return out;
}
// Comparator function for qsort
int qsortCompare(const void* a, const void* b) {
  int a_value = *static_cast<const int*>(a);
  int b_value = *static_cast<const int*>(b);
  return (a_value > b_value) - (a_value < b_value);
}
// Tests sort performance on randomly generated array of length `n_elements`
// Compares std::qsort, std::sort, quickSort and mmQuickSort algoritms
SortTestResult testRandom(int n_elements) {
  std::mt19937 mt(1337);
  Vector<int> v_qsort = test::randomIntVector(0, n_elements, n_elements,
                                                   mt);
  Vector<int> v_std_sort = v_qsort;
  Vector<int> v_quick_sort = v_qsort;
  Vector<int> v_mm_quick_sort = v_qsort;
  SortTestResult result;
  Timer timer;
  timer.start();
  std::qsort(&v_qsort[0], std::size_t(n_elements), sizeof(int), qsortCompare);
  result.qsort_time = timer.elapsed();

  timer.start();
  std::sort(v_std_sort.begin(), v_std_sort.end());
  result.std_sort_time = timer.elapsed();

  timer.start();
  quickSort(v_quick_sort.begin(), v_quick_sort.end());
  result.quick_sort_time = timer.elapsed();

  timer.start();
  mmQuickSort(v_mm_quick_sort.begin(), v_mm_quick_sort.end());
  result.mm_quick_sort_time = timer.elapsed();
  return result;
}
// Runs sort performance test `n_times` to reduce random variation 
// in the tests.
// Returns the mean of the individual results.
SortTestResult testRandomRepeat(int n_elements, int n_times) {
  SortTestResult total_time;
  total_time = total_time + testRandom(n_elements);
  return total_time / n_times;
}

}  // namespace
// Prints performance test results to std::cout
void generateReport() {
  int n = 1;
  std::cout << std::string(20, '*');
  std::cout << "Generating performance report for sorting functions" << std::endl;
  std::cout << "Sorting random vectors:" << std::endl;
  for(int i = 0; i <= 6; ++i) {
    SortTestResult result = testRandomRepeat(n, 100);
    std::cout << "n = " << n << std::endl;
    std::cout << result << std::endl;
    n *= 10;
  }
}
