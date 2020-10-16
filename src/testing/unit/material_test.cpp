#include "../../material.h"

#include <cmath>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../../geometry.h"
#include "../test_utils.h"
using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
// Test Material.apply on diffuse material when 
// light is shone directly on it
TEST(Material, ApplyOpaqueDiffuseOnPerpendicularLight) {
  Material m = Material();
  m.diffuse = Vec3(0.1, 0.2, 0.3);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector(0.0, 1.0, 0.0);
    in_vector /= in_vector.norm();
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    ASSERT_THAT(out_color, VecEq(Vec3(0.1, 0.2, 0.3) / kPi));
  }
}
TEST(Material, ApplyOpaqueDiffuseOnObliqueLight) {
  Material m = Material();
  m.diffuse = Vec3(0.1, 0.2, 0.3);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    Vec3 correct = Vec3(0.1, 0.2, 0.3) / kPi * normal.dot(in_vector);
    ASSERT_THAT(out_color, VecEq(correct));
  }
}
TEST(Material, ApplyTransparentDiffuseOnObliqueLight) {
  Material m = Material();
  m.diffuse = Vec3(0.1, 0.2, 0.3);
  m.transparent = Vec3(1.0);
  Vec3 normal(0.0, 1.0, 0.0);
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_color(1.0, 1.0, 1.0);
    Vec3 out_color = m.apply(in_vector, normal, out_vector, in_color);
    Vec3 correct = Vec3(0.1, 0.2, 0.3) / (2 * kPi) * normal.dot(in_vector);
    ASSERT_THAT(out_color, VecEq(correct));
  }
}
// Test that the sampled values never have pdf value of 0
TEST(Material, ImportanceSamplePdfPositive) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 100000; ++i) {
    Material m = test::randomMaterial(mt);

    Vec3 normal = test::randomVec3(-1.0, 1.0, mt);
    normal /= normal.norm();
    Vec3 out_vector = uniformRandomHemispherePoint(normal, mt);
    Vec3 in_vector = m.importanceSample(normal, out_vector, mt);
    double pdf = m.importanceSamplePdf(in_vector, normal, out_vector);
    ASSERT_GE(pdf, EPS);
  }
}
// Test that the importance sample pdf integrates to 1
TEST(Material, ImportanceSamplePdfIntegral) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Material m = test::randomMaterial(mt);

    Vec3 normal = test::randomVec3(0.01, 1.0, mt);
    Vec3 out_vector = test::randomVec3(0.01, 1.0, mt);
    normal /= normal.norm();
    out_vector /= out_vector.norm();
    double random_sphere_pdf = 1 / (4 * kPi);
    double sum = 0;
    int n_samples = 30000;
    for(int j = 0; j < n_samples; ++j) {
      Vec3 in_vector = uniformRandomSpherePoint(Vec3(1.0, 0.0, 0.0), mt);
      sum += m.importanceSamplePdf(in_vector, normal, out_vector);
    }
    double integral = sum / n_samples / random_sphere_pdf;
    ASSERT_NEAR(integral, 1.0, 0.1);
  }
}
// Test that the importanceSample samples follow importanceSamplePdf
TEST(Material, ImportanceSampleFollowsPdf) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 100; ++i) {
    Material m = test::randomMaterial(mt);
    Vec3 normal = uniformRandomSpherePoint(Vec3(1.0, 0.0, 0.0), mt);
    Vec3 out_vector = uniformRandomSpherePoint(Vec3(1.0, 0.0, 0.0), mt);
    Vec3 area_middle = uniformRandomSpherePoint(Vec3(1.0, 0.0, 0.0), mt);
    // minimum value for in_vector.dot(area_middle)
    double area_border = 0.8;
    double random_sphere_pdf = 1 / (4 * kPi);
    Vector<double> pdf_values;
    int n_samples = 30000;
    Vector<double> sample_values;
    for(int j = 0; j < n_samples; ++j) {
      Vec3 in_vector = uniformRandomSpherePoint(Vec3(1.0, 0.0, 0.0), mt);
      if (in_vector.dot(area_middle) > area_border) {
        double tmp = m.importanceSamplePdf(in_vector, normal, out_vector);
        pdf_values.pushBack(tmp / random_sphere_pdf);
      }
      else {
        pdf_values.pushBack(0.0);
      }
      in_vector = m.importanceSample(normal, out_vector, mt);
      if (in_vector.dot(area_middle) > area_border) {
        sample_values.pushBack(1.0);
      }
      else {
        sample_values.pushBack(0.0);
      }
    }
    double pdf_integral = test::mean(pdf_values);
    double pdf_max_error = 4 * test::standard_error_mean(pdf_values);
    double sample_fraction = test::mean(sample_values);
    double sample_max_error = 4 * test::standard_error_mean(sample_values);
    double max_error = std::max(pdf_max_error, sample_max_error);
    ASSERT_NEAR(pdf_integral, sample_fraction, max_error) << "with material:"
                                                          << std::endl
                                                          << m << std::endl;
  }
}
}  // namespace
