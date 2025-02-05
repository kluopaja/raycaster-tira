#ifndef RAYCASTER_UTILS_H
#define RAYCASTER_UTILS_H

#include <cassert>
#include <iostream>
#include <iterator>
#include <tuple>

namespace {

template <typename RandomAccessIterator>
inline RandomAccessIterator selectK(RandomAccessIterator begin,
                                    RandomAccessIterator end, std::ptrdiff_t k);
// used to sort small arrays
// has better performance for small arrays than quickSort
template <typename RandomAccessIterator>
void insertionSort(RandomAccessIterator begin, RandomAccessIterator end) {
  for (RandomAccessIterator it = begin + 1; it < end; ++it) {
    for (RandomAccessIterator it2 = it - 1; *(it2 + 1) < *it2; --it2) {
      std::swap(*it2, *(it2 + 1));
      if (it2 == begin) {
        break;
      }
    }
  }
}
// Partially sorts the elements of the array based on whether
// they are smaller, equal or greater than `*pivot`
//
// Returns pair {`n_less`, `n_same`} where `n_less` and `n_same`
// are the numbers of elements less and equal to `*pivot`
//
// Assumes that begin < end
template <typename RandomAccessIterator>
inline std::pair<std::ptrdiff_t, std::ptrdiff_t> partition(
    RandomAccessIterator begin, RandomAccessIterator end,
    RandomAccessIterator pivot) {
  std::ptrdiff_t n_same = 0;
  std::ptrdiff_t n_more = 0;
  // Every element in the interval [`begin`, `add_pos`) should be
  // in the correct place
  RandomAccessIterator add_pos = begin;
  // Counts the number of elements equal to `*pivot`
  // After this every element smaller thatn `*pivot`
  // should be in the range [`begin`, `add_pos`]
  for (RandomAccessIterator it = begin; it != end; ++it) {
    n_more += (*pivot < *it);
    if (*it < *pivot) {
      // Check if the current element is already in the correct place
      if (it == add_pos) {
        ++add_pos;
      }
      // Otherwise put it to the correct place
      else {
        // Note that now `*add_pos` will never contain a value
        // smaller than `*pivot`. If that was the case, we would
        // have already added it the range [`begin`, `add_pos`) earlier
        std::swap(*add_pos, *it);
        // Make sure that `*pivot` still contains the pivot value
        // after the swap
        if (pivot == add_pos) {
          pivot = it;
        }
        ++add_pos;
      }
    }
  }
  std::ptrdiff_t n_less = add_pos - begin;
  n_same = end - add_pos - n_more;
  if (n_same == 0) return std::make_pair(n_less, n_same);
  // now that the range [begin, begin + n_less) is filled with correct elements
  // we can no longer encounter elements for which *it < *pivot
  for (RandomAccessIterator it = add_pos + n_same; it != end;) {
    // just advance if the current element is already in the correct place
    if (*pivot < *it) {
      ++it;
    }
    // otherwise swap it to the correct place and make sure that
    // pivot stil pointing to an element equal to pivot
    else {
      std::swap(*add_pos, *it);
      pivot = add_pos++;
    }
  }
  return std::make_pair(n_less, n_same);
}
// finds median of maximum 5 elements
// uses at most 6 comparisons
template <typename RandomAccessIterator>
RandomAccessIterator medianMax5(RandomAccessIterator begin,
                                RandomAccessIterator end) {
  assert(begin < end);
  assert(end - begin <= 5);
  if (end - begin == 5) {
    // here in comments we treat the numbers as
    if (begin[1] < begin[0]) std::swap(begin[0], begin[1]);
    if (begin[3] < begin[2]) std::swap(begin[2], begin[3]);
    // 0 <= 1
    // 2 <= 3
    if (begin[2] < begin[0]) {
      std::swap(begin[0], begin[2]);
      std::swap(begin[1], begin[3]);
    }
    // possible orders:
    // 0123
    // 0213
    // 0231
    //
    // now we can ignore the 0 because it can never be the median
    // and will always be located before the median
    //
    // the median will be the third largest of 1 2 3 4
    // note that this is similar to the end - begin == 4!
    if (begin[4] < begin[1]) std::swap(begin[1], begin[4]);
    if (begin[2] < begin[1]) {
      std::swap(begin[1], begin[2]);
      std::swap(begin[3], begin[4]);
    }
    // 1 < 4, 2 < 3, 1 < 2
    // possible orders
    // 1423
    // 1243
    // 1234
    if (begin[2] < begin[4]) return begin + 2;
    return begin + 4;
  }
  if (end - begin == 4) {
    if (begin[1] < begin[0]) std::swap(begin[0], begin[1]);
    if (begin[3] < begin[2]) std::swap(begin[2], begin[3]);
    // 0 <= 1
    // 2 <= 3
    if (begin[2] < begin[0]) {
      std::swap(begin[0], begin[2]);
      std::swap(begin[1], begin[3]);
    }
    // xx = 01, yy = 23
    // xxyy
    // xyxy
    // xyyx
    if (begin[1] < begin[2]) return begin + 1;
    return begin + 2;
  }
  if (end - begin == 3) {
    if (begin[1] < begin[0]) std::swap(begin[0], begin[1]);
    // 0 < 1 2
    if (begin[2] < begin[1]) std::swap(begin[1], begin[2]);
    // 0 < 2
    // 1 < 2
    if (begin[1] < begin[0]) return begin;
    return begin + 1;
  }
  if (end - begin == 2) {
    if (begin[0] < begin[1]) return begin;
    return begin + 1;
  }
  // end - begin == 1
  return begin;
}

template <typename RandomAccessIterator>
RandomAccessIterator medianOfMedians(RandomAccessIterator begin,
                                     RandomAccessIterator end) {
  assert(begin < end);
  if (end - begin <= 5) {
    return medianMax5(begin, end);
  }
  // divide the interval into chunks of max 5
  // store the medians of these at the beginning of the container
  std::ptrdiff_t n = end - begin;
  for (std::ptrdiff_t i = 0; 5 * i < n; ++i) {
    RandomAccessIterator chunk_end = begin + n;
    if ((i + 1) * 5 < n) chunk_end = begin + (i + 1) * 5;
    RandomAccessIterator med = medianMax5(begin + i * 5, chunk_end);
    std::swap(begin[i], *med);
  }
  std::ptrdiff_t new_length = (n + 4) / 5;
  return selectK(begin, begin + new_length, (new_length - 1) / 2);
}

template <typename RandomAccessIterator>
RandomAccessIterator selectK(RandomAccessIterator begin,
                             RandomAccessIterator end, std::ptrdiff_t k) {
  assert(begin < end);
  std::ptrdiff_t n = end - begin;
  assert(k < n);
  if (n == 1) {
    return begin;
  }
  // calculate approximate median
  RandomAccessIterator median_of_medians = medianOfMedians(begin, end);
  std::ptrdiff_t n_less, n_same;
  std::tie(n_less, n_same) = partition(begin, end, median_of_medians);
  if (k >= n_less && k < n_less + n_same) {
    return begin + n_less;
  }
  if (k < n_less) {
    return selectK(begin, begin + n_less, k);
  }
  return selectK(begin + n_less + n_same, end, k - n_less - n_same);
}
// Sorts elements in range [begin, end) using quicksort with
// the median of three rule.
// Worst case O(n^2) time complexity
template <typename RandomAccessIterator>
void quickSort(RandomAccessIterator begin, RandomAccessIterator end) {
  if (end - begin <= 1) {
    return;
  }
  if (end - begin < 100) {
    insertionSort(begin, end);
    return;
  }

  // use the median of three pivot rule
  RandomAccessIterator mid = begin + (end - begin - 1) / 2;
  RandomAccessIterator pivot;
  if (*mid < *begin) std::swap(*begin, *mid);
  if (*(end - 1) < *mid) std::swap(*mid, *(end - 1));
  if (*mid < *begin)
    pivot = begin;
  else
    pivot = mid;

  std::ptrdiff_t n_less, n_same;
  std::tie(n_less, n_same) = partition(begin, end, pivot);
  quickSort(begin, begin + n_less);
  quickSort(begin + n_less + n_same, end);
}
// Sorts elements in the range [begin, end) with quicksort
// using the median of medians pivot rule.
// Worst case O(n log n) time complexity
template <typename RandomAccessIterator>
void mmQuickSort(RandomAccessIterator begin, RandomAccessIterator end) {
  if (end - begin <= 1) {
    return;
  }
  if (end - begin < 100) {
    insertionSort(begin, end);
    return;
  }

  RandomAccessIterator pivot = selectK(begin, end, (end - begin - 1) / 2);
  std::ptrdiff_t n_less, n_same;
  std::tie(n_less, n_same) = partition(begin, end, pivot);
  mmQuickSort(begin, begin + n_less);
  mmQuickSort(begin + n_less + n_same, end);
}

}  // namespace

#endif
