#ifndef RAYCASTER_VECTOR_H
#define RAYCASTER_VECTOR_H
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iterator>
#include <new>
#include <ostream>
// A class implementing a dynamic array
template <class T>
class Vector {
 public:
  using value_type = T;
  using const_iterator = const T*;
  using iterator = T*;
  // Construct a vector with 0 elements
  // Reserves memory for exactly 0 elements
  Vector() = default;
  // Constructs a vector with 'n_elements' elements of value-initialized objects
  // of type T
  //
  // Reserves memory for exactly 'n_elements' elements
  //
  // Time complexity O(n_elements)
  Vector(std::size_t n_elements);
  // Constructs a vector with elements from the initializer list `v`
  //
  // Reserves memory for exactly `v.size()` elements.
  //
  // Time complexity O(`v.size()`)
  Vector(std::initializer_list<T> v);
  // Will not reserve any more memory than needed to store the elements from 'a'
  //
  // Time complexity O(n_elements)
  Vector(const Vector& a);
  Vector(const_iterator begin, const_iterator end);
  // Leaves 'a' to the state of a valid empty vector
  //
  // Time complexity O(1)
  Vector(Vector&& a) noexcept;
  ~Vector();
  Vector& operator=(const Vector& a);
  Vector& operator=(Vector&& a) noexcept;
  // Time complexity O(1)
  T& operator[](std::size_t index);
  T& operator[](std::size_t index) const;
  // Returns a pointer to the first element
  //
  // Time complexity O(1)
  iterator begin();
  const_iterator begin() const;
  // Returns a pointer to the "element" following the last element
  //
  // Time complexity O(1)
  iterator end();
  const_iterator end() const;
  // Time complexity O(1)
  std::size_t size() const;
  // Reserves memory for total of 'new_n_reserved' elements.
  // Does not do anything if the Vector already has enough
  // reserved space for 'new_n_reserved' elements
  //
  // Time complexity O(new_n_reserved)
  void reserve(std::size_t new_n_reserved);
  // Returns the number of elements for which Vector has enough reserved
  // space
  //
  // Time complexity O(1)
  std::size_t capacity() const;
  // Shrinks reserved memory space to the minimum needed to store the
  // Vector elements
  //
  // Time complexity O(this->size())
  void shrinkToFit();
  // Removes all elements from the Vector and calls their destructors.
  // Leaves the amount of reserved space intact. Does not reallocate
  // memory.
  //
  // Time complexity O(this->size())
  void clear();
  // Time complexity amortized O(1)
  void pushBack(const T& value);
  // Time complexity amortized O(1)
  void pushBack(T&& value);
  // Removes element from the end of the Vector and calls the
  // element's destructor.
  // Undefined behaviour if the Vector is empty
  //
  // Time complexity O(1)
  void popBack();
  // If 'new_n_elements' > this->size(), adds value initialized
  // objects to the end of the vector until 'this->size()' == new_n_elements.
  // In case the new elements wouldn't fit into the reserved memory space,
  // reserves memory for exatly 'new_n_element' elements
  //
  // If 'new_n_elements' <= this->size(), removes elements (and calls their
  // destructors) until this->size() == new_n_elements.
  // Note that this does not shrink this->capacity() (i.e. does not
  // free the excess memory).
  //
  // Time complexity O(new_n_elements)
  void resize(std::size_t new_n_elements);
  // Time complexity O(1)
  void swap(Vector& a);

 private:
  const std::size_t kGrowthFactor = 2;
  std::size_t n_elements = 0;
  std::size_t n_reserved = 0;
  T* data = nullptr;
  // ensures that the vector has space for at least one
  // extra element
  void promiseOneFree();
  // As this->reserve but allows shrinking the reserved space.
  // Should not be called with new_n_reserved < this->count()!
  void forceReserve(std::size_t new_n_reserved);
  // Calls value initialization for elements in the range [start, end)
  // The memory at [start, end) should be free
  void valueInitializeRange(std::size_t start, std::size_t end);
  // Calls destructor for elements in the range [start, end)
  // The memory at [start, end) should contain valid objects
  void destructRange(std::size_t start, std::size_t end);
  // Allocates memory for n objects of type T.
  // Returns a pointer to the beginning of the allocated memory block.
  T* allocateMemory(std::size_t n);
  // Frees memory allocated by allocateMemory
  // Does nothing if 'position' == nullptr
  void freeMemory(T* position);
};
template <class T>
Vector<T>::Vector(std::size_t n_elements)
    : n_elements(n_elements), n_reserved(n_elements) {
  data = allocateMemory(n_reserved);
  valueInitializeRange(0, n_elements);
}
template <class T>
Vector<T>::Vector(std::initializer_list<T> v) {
  n_elements = v.size();
  n_reserved = v.size();
  data = allocateMemory(n_reserved);
  T* add_pos = data;
  for (auto x : v) {
    new (add_pos++) T(x);
  }
}
template <class T>
Vector<T>::Vector(const Vector& a)
    : n_elements(a.n_elements), n_reserved(a.n_elements) {
  data = allocateMemory(n_reserved);
  for (std::size_t i = 0; i < n_elements; ++i) {
    new (data + i) T(a.data[i]);
  }
}
template <class T>
Vector<T>::Vector(const_iterator begin, const_iterator end) {
  assert(std::distance(begin, end) >= 0);
  n_elements = std::size_t(std::distance(begin, end));
  n_reserved = n_elements;
  data = allocateMemory(n_reserved);
  T* add_pos = data;
  for (auto it = begin; it != end; ++it) {
    new (add_pos++) T(*it);
  }
}
template <class T>
Vector<T>::Vector(Vector&& a) noexcept
    : n_elements(a.n_elements), n_reserved(a.n_reserved), data(a.data) {
  a.n_elements = 0;
  a.n_reserved = 0;
  a.data = nullptr;
}
template <class T>
Vector<T>::~Vector() {
  destructRange(0, n_elements);
  freeMemory(data);
}
template <class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& a) {
  Vector<T> tmp(a);
  swap(tmp);
  return *this;
}
template <class T>
Vector<T>& Vector<T>::operator=(Vector&& a) noexcept {
  Vector tmp(std::move(a));
  swap(tmp);
  return *this;
}
template <class T>
inline T& Vector<T>::operator[](std::size_t index) {
  assert(index < n_elements);
  return data[index];
}
template <class T>
inline T& Vector<T>::operator[](std::size_t index) const {
  assert(index < n_elements);
  return data[index];
}
template <class T>
inline typename Vector<T>::iterator Vector<T>::begin() {
  return data;
}
template <class T>
inline typename Vector<T>::const_iterator Vector<T>::begin() const {
  return data;
}
template <class T>
inline typename Vector<T>::iterator Vector<T>::end() {
  return data + n_elements;
}
template <class T>
inline typename Vector<T>::const_iterator Vector<T>::end() const {
  return data + n_elements;
}
template <class T>
inline std::size_t Vector<T>::size() const {
  return n_elements;
}
template <class T>
void Vector<T>::reserve(std::size_t new_n_reserved) {
  if (new_n_reserved <= n_reserved) {
    return;
  }
  forceReserve(new_n_reserved);
}
template <class T>
inline std::size_t Vector<T>::capacity() const {
  return n_reserved;
}
template <class T>
void Vector<T>::shrinkToFit() {
  assert(n_elements <= n_reserved);
  if (n_elements == n_reserved) {
    return;
  }
  forceReserve(n_elements);
}
template <class T>
void Vector<T>::clear() {
  destructRange(0, n_elements);
  n_elements = 0;
}
// Currently makes a copy of 'value' to 'tmp'
// This fixes the issue where the 'value' reference gets invalidated
// during the memory reallocation step
template <class T>
inline void Vector<T>::pushBack(const T& value) {
  T tmp(value);
  promiseOneFree();
  assert(n_reserved > n_elements);
  new (data + n_elements) T(std::move(tmp));
  ++n_elements;
}
// Currently moves the 'value' to 'tmp'
// This fixes the issue where the 'value' reference gets invalidated
// during the memory reallocation step
template <class T>
inline void Vector<T>::pushBack(T&& value) {
  T tmp(std::move(value));
  promiseOneFree();
  assert(n_reserved > n_elements);
  new (data + n_elements) T(std::move(tmp));
  ++n_elements;
}
template <class T>
inline void Vector<T>::popBack() {
  assert(n_elements > 0);
  data[n_elements - 1].~T();
  --n_elements;
}
template <class T>
void Vector<T>::resize(std::size_t new_n_elements) {
  if (new_n_elements <= n_elements) {
    destructRange(new_n_elements, n_elements);
    n_elements = new_n_elements;
  } else {
    reserve(new_n_elements);
    valueInitializeRange(n_elements, new_n_elements);
    n_elements = new_n_elements;
  }
}
template <class T>
inline void Vector<T>::swap(Vector& a) {
  std::swap(n_elements, a.n_elements);
  std::swap(n_reserved, a.n_reserved);
  std::swap(data, a.data);
}
// 'n_reserved' + 1 is to handle the case 'n_reserved' == 0
template <class T>
inline void Vector<T>::promiseOneFree() {
  assert(n_reserved >= n_elements);
  if (n_elements == n_reserved) {
    reserve(kGrowthFactor * (n_reserved + 1));
  }
}
template <class T>
void Vector<T>::forceReserve(std::size_t new_n_reserved) {
  assert(new_n_reserved >= n_elements);
  T* new_data = allocateMemory(new_n_reserved);
  for (std::size_t i = 0; i < n_elements; ++i) {
    new (new_data + i) T(std::move(data[i]));
  }
  freeMemory(data);
  data = new_data;
  n_reserved = new_n_reserved;
}
template <class T>
void Vector<T>::valueInitializeRange(std::size_t start, std::size_t end) {
  for (std::size_t i = start; i < end; ++i) {
    new (data + i) T();
  }
}
template <class T>
void Vector<T>::destructRange(std::size_t start, std::size_t end) {
  for (std::size_t i = start; i < end; ++i) {
    data[i].~T();
  }
}
template <class T>
inline T* Vector<T>::allocateMemory(std::size_t n) {
  return static_cast<T*>(::operator new[](n * sizeof(T)));
}
template <class T>
inline void Vector<T>::freeMemory(T* position) {
  ::operator delete[](position);
}
template <class T>
std::ostream& operator<<(std::ostream& out, const Vector<T>& v) {
  out << "Vector(";
  for (std::size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i + 1 != v.size()) {
      out << ", ";
    }
  }
  out << ")";
  return out;
}
template <class T>
bool operator==(const Vector<T>& a, const Vector<T>& b) {
  if (a.size() != b.size()) {
    return 0;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!(a[i] == b[i])) {
      return 0;
    }
  }
  return 1;
}
#endif
