#include "raycaster.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>

using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
TEST(TestCamera, CameraInOrigo) {
  Camera camera(Vec3(0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), kPi / 2,
                kPi / 2);

  Ray r = camera.rayFromImagePlane(0.5, 0.5);
  EXPECT_THAT(r.origin, VecEq(Vec3(1.0, 0.0, 0.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, 0.0, 0.0)));
  r = camera.rayFromImagePlane(0.0, 0.0);
  EXPECT_THAT(r.origin, VecEq(Vec3(1.0, 1.0, -1.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, 1.0, -1.0)));
  r = camera.rayFromImagePlane(1.0, 1.0);
  EXPECT_THAT(r.origin, VecEq(Vec3(1.0, -1.0, 1.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, -1.0, 1.0)));
}
TEST(TestCamera, CameraNotInOrigo) {
  Camera camera(Vec3(1.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), kPi / 2,
                kPi / 2);
  Ray r = camera.rayFromImagePlane(0.5, 0.5);
  EXPECT_THAT(r.origin, VecEq(Vec3(2.0, 1.0, 1.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, 0.0, 0.0)));

  r = camera.rayFromImagePlane(0.0, 0.0);
  EXPECT_THAT(r.origin, VecEq(Vec3(2.0, 2.0, 0.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, 1.0, -1.0)));
  r = camera.rayFromImagePlane(1.0, 1.0);
  EXPECT_THAT(r.origin, VecEq(Vec3(2.0, 0.0, 2.0)));
  EXPECT_THAT(r.direction, VecEq(Vec3(1.0, -1.0, 1.0)));
}

}  // namespace
