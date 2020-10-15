#include "../../lighting_utils.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../test_utils.h"

using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
TEST(PerfectRefraction, NormalVector) {
  Vec3 in_vector(0.0, 1.0, 0.0);
  Vec3 normal(0.0, 1.0, 0.0);
  double eta_1 = 1.1;
  double eta_2 = 2.3;
  Vec3 out_vector;
  bool success;
  std::tie(out_vector, success) = perfectRefraction(in_vector, normal,
                                                    eta_1, eta_2);
  Vec3 correct(0.0, -1.0, 0.0);
  EXPECT_THAT(out_vector, VecEq(correct));
  EXPECT_TRUE(success);
}
TEST(PerfectRefraction, SameRefractiveIndexSimple) {
  Vec3 in_vector(-1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0), 0.0);
  Vec3 normal(0.0, 1.0, 0.0);
  double eta_1 = 4.0;
  double eta_2 = 4.0;
  Vec3 out_vector;
  bool success;
  std::tie(out_vector, success) = perfectRefraction(in_vector, normal,
                                                    eta_1, eta_2);
  Vec3 correct(1.0/std::sqrt(2.0), -1.0/std::sqrt(2.0), 0.0);
  EXPECT_THAT(out_vector, VecEq(correct));
  EXPECT_TRUE(success);
}
TEST(PerfectRefraction, SameRefractiveIndex) {
  Vec3 in_vector(1.0, 4.0, 0.0);
  Vec3 normal(1.1, 2.2, 3.3);
  double eta_1 = 4.0;
  double eta_2 = 4.0;
  Vec3 out_vector;
  bool success;
  std::tie(out_vector, success) = perfectRefraction(in_vector, normal,
                                                    eta_1, eta_2);
  Vec3 correct(-1.0, -4.0, 0.0);
  EXPECT_THAT(out_vector, VecEq(correct));
  EXPECT_TRUE(success);
}
TEST(PerfectRefraction, SameRefractiveIndexRandom) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 10000; ++i) {
    Vec3 in_vector = test::randomVec3(-2.0, 2.0, mt);
    Vec3 normal = test::randomVec3(-2.0, 2.0, mt);
    double eta_1 = 4.0;
    double eta_2 = 4.0;
    Vec3 out_vector;
    bool success;
    std::tie(out_vector, success) = perfectRefraction(in_vector, normal,
        eta_1, eta_2);
    ASSERT_TRUE(success);
    ASSERT_THAT(out_vector, VecEq(-1.0 * in_vector));
    Vec3 out_vector_2;
    std::tie(out_vector_2, success) = perfectRefraction(out_vector, normal,
        eta_1, eta_2);
    ASSERT_TRUE(success);
    ASSERT_THAT(in_vector, VecEq(out_vector_2));
  }
}
// test that refraction is reversible if it is successful
TEST(PerfectRefraction, IsReversibleRandom) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 10000; ++i) {
    Vec3 in_vector = test::randomVec3(-2.0, 2.0, mt);
    Vec3 normal = test::randomVec3(-2.0, 2.0, mt);
    std::uniform_real_distribution dist(0.1, 10.0);
    double eta_1 = dist(mt);
    double eta_2 = dist(mt);
    Vec3 out_vector;
    bool success;
    std::tie(out_vector, success) = perfectRefraction(in_vector, normal,
                                                      eta_1, eta_2);
    if(!success) continue;
    Vec3 out_vector_2;
    std::tie(out_vector_2, success) = perfectRefraction(out_vector, normal,
        eta_1, eta_2);
    ASSERT_THAT(in_vector, VecEq(out_vector_2));
  }
}
}
