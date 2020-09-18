#include "utils.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <chrono>
#include <random>

#include "test_utils.h"

using ::testing::ContainerEq;
using ::testing::PrintToString;
using ::testing::UnorderedElementsAreArray;

namespace {

TEST(InsertionSort, Simple) {
  std::vector<int> v = {5, 4, 3, 2, 1};
  insertionSort(v.begin(), v.end());
  std::vector<int> c = {1, 2, 3, 4, 5};
  EXPECT_THAT(v, ContainerEq(c));
}
TEST(InsertionSort, Random) {
  std::mt19937 mt(1337);
  int n_tests = 100'000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(1, 30);
    int n_values = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 30, n_values, mt);
    std::vector<int> v_copy = v;
    std::sort(v_copy.begin(), v_copy.end());
    insertionSort(v.begin(), v.end());
    ASSERT_THAT(v, ContainerEq(v_copy));
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(MedianMax5, Random) {
  std::mt19937 mt(1337);
  int n_tests = 100'000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution length_dist(1, 5);
    std::uniform_int_distribution max_value_dist(1, 30);
    int max_value = max_value_dist(mt);
    int n_values = length_dist(mt);
    std::vector<int> v = test::randomIntVector(0, max_value, n_values, mt);
    int median = *medianMax5(v.begin(), v.end());
    std::sort(v.begin(), v.end());
    ASSERT_EQ(median, v[(n_values - 1) / 2])
        << " of values " + PrintToString(v) << std::endl;
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(MedianOfMedians, Random) {}
TEST(CountRelativeValues, Simple) {
  std::vector<int> v = {0, 1, 2, 3, 3, 4, 5};
  size_t n_less, n_same;
  std::tie(n_less, n_same) =
      countRelativeValues(v.begin(), v.end(), v.begin() + 3);
  EXPECT_THAT(n_less, 3);
  EXPECT_THAT(n_same, 2);
}
TEST(Partition, Simple) {
  std::vector<int> v = {3, 2, 2, 1, 1, 0};
  partition(v.begin(), v.end(), 3, 2, v.begin() + 1);
  std::vector<int> p1(v.begin(), v.begin() + 3);
  std::vector<int> p2(v.begin() + 3, v.begin() + 5);
  std::vector<int> p3(v.begin() + 5, v.end());
  std::vector<int> c1 = {0, 1, 1};
  std::vector<int> c2 = {2, 2};
  std::vector<int> c3 = {3};
  EXPECT_THAT(p1, UnorderedElementsAreArray(c1.begin(), c1.end()));
  EXPECT_THAT(p2, UnorderedElementsAreArray(c2.begin(), c2.end()));
  EXPECT_THAT(p3, UnorderedElementsAreArray(c3.begin(), c3.end()));
}
TEST(Partition, Random) {
  std::mt19937 mt(1337);
  int n_tests = 1000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(1, 100);
    int n = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 100, n, mt);
    std::uniform_int_distribution pivot_dist(0, n - 1);
    std::vector<int>::iterator pivot = v.begin() + pivot_dist(mt);

    int n_less, n_same;
    std::tie(n_less, n_same) = countRelativeValues(v.begin(), v.end(), pivot);
    // calculate the correct result
    std::vector<int> v_copy = v;
    std::sort(v_copy.begin(), v_copy.end());
    std::vector<int> c1(v_copy.begin(), v_copy.begin() + n_less);
    std::vector<int> c2(v_copy.begin() + n_less,
                        v_copy.begin() + n_less + n_same);
    std::vector<int> c3(v_copy.begin() + n_less + n_same, v_copy.end());

    partition(v.begin(), v.end(), n_less, n_same, pivot);
    std::vector<int> p1(v.begin(), v.begin() + n_less);
    std::vector<int> p2(v.begin() + n_less, v.begin() + n_less + n_same);
    std::vector<int> p3(v.begin() + n_less + n_same, v.end());

    ASSERT_THAT(p1, UnorderedElementsAreArray(c1.begin(), c1.end()));
    ASSERT_THAT(p2, UnorderedElementsAreArray(c2.begin(), c2.end()));
    ASSERT_THAT(p3, UnorderedElementsAreArray(c3.begin(), c3.end()));
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(SelectK, Simple) {
  std::vector<int> v = {0, 1, 2, 3};
  int kth_value = *selectK(v.begin(), v.end(), 0);
  EXPECT_EQ(kth_value, 0) << " of values " << PrintToString(v) << std::endl;

  v = {0, 0, 0, 0};
  kth_value = *selectK(v.begin(), v.end(), 3);
  EXPECT_EQ(kth_value, 0) << " of values " << PrintToString(v) << std::endl;

  v = {4, 3, 2, 1};
  kth_value = *selectK(v.begin(), v.end(), 2);
  EXPECT_EQ(kth_value, 3) << " of values " << PrintToString(v) << std::endl;
}
TEST(SelectK, Random) {
  std::mt19937 mt(1337);
  int n_tests = 100'000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution length_dist(1, 100);
    std::uniform_int_distribution max_value_dist(1, 100);
    int max_value = max_value_dist(mt);
    int n_values = length_dist(mt);
    std::vector<int> v = test::randomIntVector(0, max_value, n_values, mt);
    std::uniform_int_distribution k_dist(0, n_values - 1);
    int k = k_dist(mt);
    int kth_value = *selectK(v.begin(), v.end(), k);
    std::sort(v.begin(), v.end());
    ASSERT_EQ(kth_value, v[k])
        << " of values " + PrintToString(v) << std::endl
        << "after passing " << i << " random tests" << std::endl;
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(QuickSort, Simple) {
  std::vector<int> v = {4, 3, 2, 1};
  quickSort(v.begin(), v.end());
  std::vector<int> c = {1, 2, 3, 4};
  EXPECT_THAT(v, ContainerEq(c));

  v = {0, 2, 1, 3};
  quickSort(v.begin(), v.end());
  c = {0, 1, 2, 3};
  EXPECT_THAT(v, ContainerEq(c));
}
TEST(QuickSort, RandomSmall) {
  std::mt19937 mt(1337);
  int n_tests = 100'000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(1, 30);
    int n_values = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 30, n_values, mt);
    std::vector<int> v_copy = v;
    quickSort(v.begin(), v.end());
    std::sort(v_copy.begin(), v_copy.end());
    ASSERT_THAT(v, ContainerEq(v_copy));
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(QuickSort, RandomLarge) {
  std::mt19937 mt(1337);
  int n_tests = 50;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(1, 1'000'000);
    int n_values = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 1'000'000, n_values, mt);
    std::vector<int> v_copy = v;
    quickSort(v.begin(), v.end());
    std::sort(v_copy.begin(), v_copy.end());
    ASSERT_THAT(v, ContainerEq(v_copy));
  }
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(QuickSort, LargePerformance) {
  std::mt19937 mt(1337);
  std::chrono::duration<double, std::milli> std_sort_time(0.0);
  std::chrono::duration<double, std::milli> quick_sort_time(0.0);
  int n_tests = 50;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(100'000, 1'000'000);
    int n_values = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 1'000'000, n_values, mt);
    std::vector<int> v_copy = v;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::sort(v_copy.begin(), v_copy.end());
    auto t1 = std::chrono::high_resolution_clock::now();
    quickSort(v.begin(), v.end());
    auto t2 = std::chrono::high_resolution_clock::now();

    std_sort_time += t1 - t0;
    quick_sort_time += t2 - t1;
    ASSERT_THAT(v, ContainerEq(v_copy));
  }
  EXPECT_LE(quick_sort_time.count(), 3 * std_sort_time.count())
      << "Did not reach the performance requirements";
  std::cerr << "std_sort_time: " << std_sort_time.count() << std::endl;
  std::cerr << "quick_sort_time: " << quick_sort_time.count() << std::endl;
  std::cerr << n_tests << " random test cases run" << std::endl;
}
TEST(QuickSort, SmallPerformance) {
  std::mt19937 mt(1337);
  std::chrono::duration<double, std::milli> std_sort_time(0.0);
  std::chrono::duration<double, std::milli> quick_sort_time(0.0);
  int n_tests = 100'000;
  for (int i = 0; i < n_tests; ++i) {
    std::uniform_int_distribution dist(1, 200);
    int n_values = dist(mt);
    std::vector<int> v = test::randomIntVector(1, 200, n_values, mt);
    std::vector<int> v_copy = v;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::sort(v_copy.begin(), v_copy.end());
    auto t1 = std::chrono::high_resolution_clock::now();
    quickSort(v.begin(), v.end());
    auto t2 = std::chrono::high_resolution_clock::now();

    std_sort_time += t1 - t0;
    quick_sort_time += t2 - t1;
    ASSERT_THAT(v, ContainerEq(v_copy));
  }
  EXPECT_LE(quick_sort_time.count(), 3 * std_sort_time.count())
      << "Did not reach the performance requirements";
  std::cerr << "std_sort_time: " << std_sort_time.count() << std::endl;
  std::cerr << "quick_sort_time: " << quick_sort_time.count() << std::endl;
  std::cerr << n_tests << " random test cases run" << std::endl;
}

}  // namespace
