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
MATCHER_P2(VecNear, v, e, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < e;
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
TEST(RenderTest, EnvironmentLightColor) {
  Camera camera(Vec3(0, 0.0, .0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 5.0, kPi / 5.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  Image result = scene.render(50, 50, 10, 3, 1);
  EXPECT_THAT(result.getColor(10, 0), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(10, 40), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(25, 25), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(49, 49), VecEq(Vec3(1.0)));
}
// The furnace test
// https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
// Every point of the ball is illuminated equally from every direction
// (from the outer side) with intensity 1.
// The ball will reflect 0.18 of the received light to every direction
// with equal radiance.
// Therefore the ball should look like a disk with 0.18 intensity
//
// Furnace test with uniform sphere sampling
TEST(RenderTest, OpaqueDiffuseUniformSampling) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 7.0, kPi / 7.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.addModelFromFile("../models/testing/ball_diffuse.obj", Vec3(0.0), kSmooth);
  scene.setSamplingScheme(kUniformSphere);
  Image result = scene.render(20, 20, 1000, 4, 1);
  result.savePPM("test.ppm");
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(0.18), 0.03));
  EXPECT_THAT(result.getColor(0, 19), VecNear(Vec3(0.18), 0.03));
  EXPECT_THAT(result.getColor(10, 10), VecNear(Vec3(0.18), 0.03));
  EXPECT_THAT(result.getColor(9, 9), VecNear(Vec3(0.18), 0.03));
}
// Furnace test with importance sampling
TEST(RenderTest, OpaqueDiffuseImportanceSampling) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 7.0, kPi / 7.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.addModelFromFile("../models/testing/ball_diffuse.obj", Vec3(0.0), kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(20, 20, 100, 4, 1);
  result.savePPM("test.ppm");
  // With importance sampling the result should be very close to 0.18
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(0, 19), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(10, 10), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(9, 9), VecNear(Vec3(0.18), 0.001));
}
// Furnace test with transparent diffuse
TEST(RenderTest, TransparentDiffuseUniformSampling) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 7.0, kPi / 7.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.addModelFromFile("../models/testing/ball_transparent_diffuse.obj", Vec3(0.0), kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(20, 20, 100, 4, 1);
  result.savePPM("test.ppm");
  // With importance sampling the result should be very close to 0.18
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(0, 19), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(10, 10), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(9, 9), VecNear(Vec3(0.18), 0.001));
}
}  // namespace
