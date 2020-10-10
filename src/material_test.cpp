#include "material.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.h"
#include "test_utils.h"
using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
// Test Material.apply on diffuse material when 
// light is shone directly on it
TEST(Material, ApplyOpaqueDiffuseOnPerpendicularLight) {
  Material m = Material();
  m.diffuse = Vec3(1.0, 1.0, 1.0);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector(0.0, 1.0, 0.0);
    in_vector /= in_vector.norm();
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    ASSERT_THAT(out_color, VecEq(Vec3(1.0) / kPi));
  }
}
TEST(Material, ApplyOpaqueDiffuseOnObliqueLight) {
  Material m = Material();
  m.diffuse = Vec3(1.0, 1.0, 1.0);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    Vec3 correct = Vec3(1.0) / kPi * normal.dot(in_vector);
    ASSERT_THAT(out_color, VecEq(correct));
  }
}
TEST(Material, ApplyTransparentDiffuseOnObliqueLight) {
  Material m = Material();
  m.diffuse = Vec3(1.0, 1.0, 1.0);
  m.transparent = Vec3(1.0);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    Vec3 correct = Vec3(1.0) / (2 * kPi) * normal.dot(in_vector);
    ASSERT_THAT(out_color, VecEq(correct));
  }
}
// Test that the sampled values never have pdf value of 0
TEST(Material, ImportanceSamplePdfPositive) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 100000; ++i) {
    Material m = Material();
    m.diffuse = test::randomVec3(0.01, 1.0, mt);
    m.emitted = test::randomVec3(0.01, 1.0, mt);
    m.specular = test::randomVec3(0.01, 1.0, mt);
    std::uniform_real_distribution dist(0.1, 100.0);
    m.specular_exp = dist(mt);
    m.index_of_refraction = dist(mt);
    Vec3 normal = test::randomVec3(-1.0, 1.0, mt);
    normal /= normal.norm();
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = m.importanceSample(normal, out_vector, mt);
    double pdf = m.importanceSamplePdf(in_vector, normal, out_vector);
    ASSERT_GE(pdf, EPS);
  }
}
}  // namespace
