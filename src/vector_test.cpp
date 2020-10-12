#include "vector.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <random>
#include <sstream>
#include <string>
#include <vector>

using ::testing::PrintToString;

namespace {
MATCHER_P(ElementType1Eq, et,
          "should equal {" + PrintToString(et.a) + ", " + PrintToString(et.b)) {
  return (et.a == arg.a) && (et.b == arg.b);
}
// Class to be used as a test element in Vector
class ElementType1 {
 public:
  int a;
  int b = 2;
};
TEST(Vector, DefaultConstructor) {
  Vector<int> v;
  EXPECT_THAT(v.size(), 0u);

  Vector<ElementType1> v2;
  EXPECT_THAT(v.size(), 0u);
}
// Tests that Vector(std::size_t n_elements) creates vectors of correct
// size
TEST(Vector, SizeConstructorSize) {
  Vector<int> v(0);
  EXPECT_THAT(v.size(), 0u);

  Vector<int> v2(10);
  EXPECT_THAT(v2.size(), 10u);

  Vector<ElementType1> v3(1000);
  EXPECT_THAT(v3.size(), 1000u);
}
// TODO test the destructor
TEST(Vector, OperatorSquareBrackets) {
  Vector<int> v(10);
  v[0] = 0;
  v[1] = 1;
  v[2] = 2;
  v[9] = 9;
  EXPECT_EQ(v[0], 0);
  EXPECT_EQ(v[1], 1);
  EXPECT_EQ(v[2], 2);
  EXPECT_EQ(v[9], 9);
}
// Tests the operator[](std::size_t index) const
TEST(Vector, OperatorSquareBracketsConst) {
  Vector<int> v(10);
  v[0] = 0;
  v[1] = 1;
  v[2] = 2;
  v[9] = 9;
  const Vector<int> const_v(v);
  EXPECT_EQ(const_v[0], 0);
  EXPECT_EQ(const_v[1], 1);
  EXPECT_EQ(const_v[2], 2);
  EXPECT_EQ(const_v[9], 9);
}
// Tests that Vector(std::size_t n_elements) value initializes the
// elements correctly
TEST(Vector, SizeConstructorInitialization) {
  Vector<int> v1(10);
  ASSERT_EQ(v1.size(), 10u);
  for (std::size_t i = 0; i < 10; ++i) {
    ASSERT_EQ(v1[i], 0);
  }
  Vector<ElementType1> v2(100);
  ASSERT_THAT(v2.size(), 100u);
  for (std::size_t i = 0; i < 100; ++i) {
    ASSERT_THAT(v2[i].a, 0);
    ASSERT_THAT(v2[i].b, 2);
  }
}
TEST(Vector, CopyConstructor) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  Vector<int> v2(v1);
  ASSERT_EQ(v2.size(), 3);
  EXPECT_EQ(v2[0], 0);
  EXPECT_EQ(v2[1], 1);
  EXPECT_EQ(v2[2], 2);
  ASSERT_EQ(v1.size(), 3);
  EXPECT_EQ(v1[0], 0);
  EXPECT_EQ(v1[1], 1);
  EXPECT_EQ(v1[2], 2);
}
// Tests that the modifications of the copy and original do not
// affect each other
TEST(Vector, CopyConstructorModify) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  Vector<int> v2(v1);
  v1[0] = 10;
  v2[0] = 20;
  EXPECT_EQ(v1[0], 10);
  EXPECT_EQ(v2[0], 20);
}
TEST(Vector, MoveConstructor) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  Vector<int> v2(std::move(v1));
  ASSERT_EQ(v2.size(), 3);
  EXPECT_EQ(v2[0], 0);
  EXPECT_EQ(v2[1], 1);
  EXPECT_EQ(v2[2], 2);
  ASSERT_EQ(v1.size(), 0);
}
TEST(Vector, CopyAssignment) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  Vector<int> v2(3);
  v2[0] = 2;
  v2[1] = 1;
  v2[2] = 0;
  v1 = v2;
  ASSERT_EQ(v1.size(), 3);
  EXPECT_EQ(v2[0], 2);
  EXPECT_EQ(v2[1], 1);
  EXPECT_EQ(v2[2], 0);
}
TEST(Vector, CopyAssignmentSelf) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  v1 = v1;
  ASSERT_EQ(v1.size(), 3);
  EXPECT_EQ(v1[0], 0);
  EXPECT_EQ(v1[1], 1);
  EXPECT_EQ(v1[2], 2);
}
TEST(Vector, MoveAssignment) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  Vector<int> v2(3);
  v2[0] = 2;
  v2[1] = 1;
  v2[2] = 0;
  v1 = std::move(v2);
  ASSERT_EQ(v2.size(), 0);
  ASSERT_EQ(v1.size(), 3);
  EXPECT_EQ(v1[0], 2);
  EXPECT_EQ(v1[1], 1);
  EXPECT_EQ(v1[2], 0);
}
TEST(Vector, Begin) {
  Vector<int> v(3);
  v[0] = 10;
  EXPECT_EQ(v[0], *v.begin());
  // Test that the iterator can be used to edit the correct element
  (*v.begin()) = 5;
  EXPECT_EQ(v[0], 5);
}
TEST(Vector, End) {
  Vector<int> v(3);
  v[2] = 10;
  int* begin = v.begin();
  int* end = v.end();
  EXPECT_EQ(std::distance(begin, end), 3);
  EXPECT_EQ(*(--end), 10);
}
TEST(Vector, Capacity) {
  Vector<int> v;
  EXPECT_EQ(v.capacity(), 0u);
  v = Vector<int>(20);
  EXPECT_EQ(v.capacity(), 20u);
}
TEST(Vector, Reserve) {
  Vector<int> v(2);
  v[0] = 1;
  v[1] = 2;
  ASSERT_EQ(v.capacity(), 2u);
  v.reserve(100);
  // Test that the elements are still correct
  ASSERT_EQ(v.size(), 2);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v.capacity(), 100u);
  // Test that the reserve doesn't shrink the reserved space
  v.reserve(0);
  EXPECT_EQ(v.capacity(), 100u);
}
TEST(Vector, ShrinkToFit) {
  Vector<int> v(2);
  v[0] = 1;
  v[1] = 2;
  v.reserve(100);
  v.shrinkToFit();
  EXPECT_EQ(v.capacity(), 2u);
  // Test that the elements are still correct
  ASSERT_EQ(v.size(), 2);
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
}
TEST(Vector, Clear) {
  Vector<int> v(2);
  v[0] = 1;
  v[1] = 2;
  v.clear();
  EXPECT_EQ(v.size(), 0);
  EXPECT_EQ(v.capacity(), 2);
}
TEST(Vector, PushBackCopy) {
  Vector<ElementType1> v;
  v.pushBack({1, 2});
  v.pushBack({2, 3});
  v.pushBack({3, 4});
  ASSERT_EQ(v.size(), 3);
  EXPECT_THAT(v[0], ElementType1Eq(ElementType1{1, 2}));
  EXPECT_THAT(v[1], ElementType1Eq(ElementType1{2, 3}));
  EXPECT_THAT(v[2], ElementType1Eq(ElementType1{3, 4}));
}
TEST(Vector, PushBackMove) {
  Vector<ElementType1> v;
  ElementType1 e1 = {1, 2};
  ElementType1 e2 = {2, 3};
  ElementType1 e3 = {3, 4};
  v.pushBack(std::move(e1));
  v.pushBack(std::move(e2));
  v.pushBack(std::move(e3));
  ASSERT_EQ(v.size(), 3);
  EXPECT_THAT(v[0], ElementType1Eq(ElementType1{1, 2}));
  EXPECT_THAT(v[1], ElementType1Eq(ElementType1{2, 3}));
  EXPECT_THAT(v[2], ElementType1Eq(ElementType1{3, 4}));
}
TEST(Vector, PopBack) {
  Vector<int> v(2);
  v[0] = 1;
  v[1] = 2;
  v.popBack();
  ASSERT_THAT(v.size(), 1u);
  EXPECT_EQ(v[0], 1);
  v.popBack();
  ASSERT_EQ(v.size(), 0u);
}
// Test increasing the Vector size with Vector.resize
TEST(Vector, ResizeIncrease) {
  Vector<int> v(2);
  v[0] = 1;
  v[1] = 2;
  v.resize(4);
  ASSERT_THAT(v.capacity(), 4u);
  ASSERT_THAT(v.size(), 4u);
  EXPECT_THAT(v[0], 1);
  EXPECT_THAT(v[1], 2);
  EXPECT_THAT(v[2], 0);
  EXPECT_THAT(v[3], 0);
}
// Test decreaseing the Vector size with Vector.resize
TEST(Vector, ResizeDecrease) {
  Vector<int> v(3);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  v.resize(1);
  ASSERT_THAT(v.capacity(), 3u);
  ASSERT_THAT(v.size(), 1u);
  EXPECT_THAT(v[0], 1);
}
TEST(Vector, Swap) {
  Vector<int> v1(3);
  v1[0] = 0;
  v1[1] = 1;
  v1[2] = 2;
  v1.reserve(10);
  Vector<int> v2(2);
  v2[0] = 3;
  v2[1] = 4;
  v1.swap(v2);
  EXPECT_THAT(v1.capacity(), 2);
  ASSERT_THAT(v1.size(), 2);
  EXPECT_THAT(v1[0], 3);
  EXPECT_THAT(v1[1], 4);
  EXPECT_THAT(v2.capacity(), 10);
  ASSERT_THAT(v2.size(), 3);
  EXPECT_THAT(v2[0], 0);
  EXPECT_THAT(v2[1], 1);
  EXPECT_THAT(v2[2], 2);
}
TEST(Vector, Random) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 1000; ++i) {
    Vector<int> v_own;
    std::vector<int> v_std;
    std::uniform_real_distribution operation_dist(0.0, 1.0);
    std::uniform_int_distribution dist_0_1000(0, 1000);
    for (int j = 0; j < 10000; ++j) {
      double operation = operation_dist(mt);
      if (operation < 0.3) {
        ASSERT_EQ(v_own.size(), v_std.size());
        if (v_std.size() > 0) {
          std::uniform_int_distribution<std::size_t> pos_dist(0,
                                                              v_std.size() - 1);
          std::size_t pos = pos_dist(mt);
          int value = dist_0_1000(mt);
          v_own[pos] = value;
          v_std[pos] = value;
        }
      } else if (operation < 0.35) {
        std::size_t tmp = std::size_t(dist_0_1000(mt));
        v_own.reserve(tmp);
        v_std.reserve(tmp);
      } else if (operation < 0.40) {
        v_own.shrinkToFit();
        v_std.shrink_to_fit();
      } else if (operation < 0.45) {
        v_own.clear();
        v_std.clear();
      } else if (operation < 0.70) {
        ASSERT_EQ(v_own.size(), v_std.size());
        if (v_std.size() > 0) {
          v_own.popBack();
          v_std.pop_back();
        }
      } else if (operation < 0.95) {
        int tmp = std::size_t(dist_0_1000(mt));
        v_own.pushBack(tmp);
        v_std.push_back(tmp);
      } else {
        std::size_t tmp = std::size_t(dist_0_1000(mt));
        v_own.resize(tmp);
        v_std.resize(tmp);
      }
    }
    ASSERT_EQ(v_own.size(), v_std.size());
    for (std::size_t i = 0; i < v_own.size(); ++i) {
      ASSERT_EQ(v_own[i], v_std[i]);
    }
  }
}
// Test operator <<
TEST(Vector, ToOstream) {
  Vector<int> v(3);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  std::stringstream ss;
  ss << v;
  std::string correct = "Vector(1, 2, 3)";
  ASSERT_TRUE(ss.str() == correct);
}
}  // namespace
