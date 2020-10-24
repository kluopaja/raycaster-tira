#include "../../raycaster.h"

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
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(1.0, 0.0, 0.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, 0.0, 0.0)));
  r = camera.rayFromImagePlane(0.0, 0.0);
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(1.0, 1.0, -1.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, 1.0, -1.0)));
  r = camera.rayFromImagePlane(1.0, 1.0);
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(1.0, -1.0, 1.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, -1.0, 1.0)));
}
TEST(TestCamera, CameraNotInOrigo) {
  Camera camera(Vec3(1.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), kPi / 2,
                kPi / 2);
  Ray r = camera.rayFromImagePlane(0.5, 0.5);
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(2.0, 1.0, 1.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, 0.0, 0.0)));

  r = camera.rayFromImagePlane(0.0, 0.0);
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(2.0, 2.0, 0.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, 1.0, -1.0)));
  r = camera.rayFromImagePlane(1.0, 1.0);
  EXPECT_THAT(r.getOrigin(), VecEq(Vec3(2.0, 0.0, 2.0)));
  EXPECT_THAT(r.getDirection(), VecEq(Vec3(1.0, -1.0, 1.0)));
}
TEST(EnvironmentLight, Constructor) {
  EnvironmentLight e;
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0)), VecEq(Vec3(0.0)));
}
TEST(EnvironmentLight, UniformLight) {
  EnvironmentLight e;
  e.setColor(Vec3(1.0, 2.0, 3.0));
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0)), VecEq(Vec3(1.0, 2.0, 3.0)));
  EXPECT_THAT(e.colorAtDirection(Vec3(0.0, 1.0, 0.0)), VecEq(Vec3(1.0, 2.0, 3.0)));
  EXPECT_THAT(e.colorAtDirection(Vec3(-1.0)), VecEq(Vec3(1.0, 2.0, 3.0)));
}
TEST(EnvironmentLight, DirectedLightOneExponent) {
  EnvironmentLight e;
  e.setColor(Vec3(1.0, 2.0, 3.));
  e.setDirected(Vec3(1.0, 0.0, 0.0), 1.0);
  // opposite direction
  EXPECT_THAT(e.colorAtDirection(Vec3(-1.0, 0.0, 0.0)), VecEq(Vec3(0.0, 0.0, 0.0)));
  // perpendicular direction should also be 0
  EXPECT_THAT(e.colorAtDirection(Vec3(0.0, 1.0, 0.0)), VecEq(Vec3(0.0, 0.0, 0.0)));
  // calculate the cos attenuation
  Vec3 correct = std::cos(kPi/4) * Vec3(1.0, 2.0, 3.);
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0, 1.0, 0.0)), VecEq(correct));
  // direct should also be 1
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0, 0.0, 0.0)), VecEq(Vec3(1.0, 2.0, 3.0)));
}
TEST(EnvironmentLight, DirectedLightLargeExponent) {
  EnvironmentLight e;
  e.setColor(Vec3(1.0, 2.0, 3.));
  e.setDirected(Vec3(1.0, 0.0, 0.0), 123.0);
  // opposite direction
  EXPECT_THAT(e.colorAtDirection(Vec3(-1.0, 0.0, 0.0)), VecEq(Vec3(0.0, 0.0, 0.0)));
  // perpendicular direction should also be 0
  EXPECT_THAT(e.colorAtDirection(Vec3(0.0, 1.0, 0.0)), VecEq(Vec3(0.0, 0.0, 0.0)));
  // calculate the cos attenuation
  Vec3 correct = std::pow(std::cos(kPi/4), 123.0) * Vec3(1.0, 2.0, 3.);
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0, 1.0, 0.0)), VecEq(correct));
  // direct should also be 1
  EXPECT_THAT(e.colorAtDirection(Vec3(1.0, 0.0, 0.0)), VecEq(Vec3(1.0, 2.0, 3.0)));
}
TEST(RenderTest, EnvironmentLightColor) {
  Camera camera(Vec3(0, 0.0, 0.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 5.0, kPi / 5.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  Image result = scene.render(50, 50, 10, 3, 1);
  EXPECT_THAT(result.getColor(10, 0), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(10, 40), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(25, 25), VecEq(Vec3(1.0)));
  EXPECT_THAT(result.getColor(49, 49), VecEq(Vec3(1.0)));
}
// Tests directed environment light with the exponent 1.0
// Renders very narrow view to various directions
TEST(RenderTest, EnvironmentLightDirectedExponent1) {

  // Opposite way from where the light is coming
  Camera camera(Vec3(0, 0.0, 0.0), Vec3(0.0, -1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 1.0);
  Image result = scene.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result.getColor(0, 0), VecEq(Vec3(0.0)));

  // Directly towards the light
  Camera camera2(Vec3(0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene2(camera2);
  scene2.setEnvironmentLightColor(Vec3(1.0));
  scene2.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 1.0);
  Image result2 = scene2.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result2.getColor(0, 0), VecEq(Vec3(1.0)));


  // 45 degress away from directly towards the light
  Camera camera3(Vec3(0, 0.0, 0.0), Vec3(1.0, 1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene3(camera3);
  scene3.setEnvironmentLightColor(Vec3(1.0));
  scene3.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 1.0);
  Image result3 = scene3.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result3.getColor(0, 0), VecEq(std::cos(kPi / 4.0) * Vec3(1.0)));
}
// Tests directed environment light with exponent5
// Renders very narrow view to various directions
TEST(RenderTest, EnvironmentLightDirectedExponent5) {

  // Opposite way from where the light is coming
  Camera camera(Vec3(0, 0.0, 0.0), Vec3(0.0, -1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 5.0);
  Image result = scene.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result.getColor(0, 0), VecEq(Vec3(0.0)));

  // Directly towards the light
  Camera camera2(Vec3(0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene2(camera2);
  scene2.setEnvironmentLightColor(Vec3(1.0));
  scene2.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 5.0);
  Image result2 = scene2.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result2.getColor(0, 0), VecEq(Vec3(1.0)));


  // 45 degress away from directly towards the light
  Camera camera3(Vec3(0, 0.0, 0.0), Vec3(1.0, 1.0, 0.0),
                Vec3(0.0, 0.0, 1.0), kPi / 100.0, kPi / 100.0);
  Scene scene3(camera3);
  scene3.setEnvironmentLightColor(Vec3(1.0));
  scene3.setEnvironmentLightDirected(Vec3(0.0, 1.0, 0.0), 5.0);
  Image result3 = scene3.render(1, 1, 1, 1, 1);
  EXPECT_THAT(result3.getColor(0, 0),
              VecEq(std::pow(std::cos(kPi / 4.0), 5.0) * Vec3(1.0)));
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
  scene.addModelFromFile("../models/testing/ball_diffuse.obj", Vec3(0.0), NormalType::kSmooth);
  scene.setSamplingScheme(kUniformSphere);
  Image result = scene.render(20, 20, 1000, 4, 1);
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
  scene.addModelFromFile("../models/testing/ball_diffuse.obj", Vec3(0.0), NormalType::kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(20, 20, 100, 4, 1);
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
  scene.addModelFromFile("../models/testing/ball_transparent_diffuse.obj", Vec3(0.0), NormalType::kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(20, 20, 100, 4, 1);
  // With importance sampling the result should be very close to 0.18
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(0, 19), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(10, 10), VecNear(Vec3(0.18), 0.001));
  EXPECT_THAT(result.getColor(9, 9), VecNear(Vec3(0.18), 0.001));
}
// Test mirror in white environment
TEST(RenderTest, MirrorEnvironment) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 1.1, kPi / 1.1);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.addModelFromFile("../models/testing/XY_mirror.obj", Vec3(0.0), NormalType::kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(100, 100, 100, 4, 1);
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(1.0), 0.001));
  EXPECT_THAT(result.getColor(99, 99), VecNear(Vec3(1.0), 0.001));
  EXPECT_THAT(result.getColor(20, 36), VecNear(Vec3(1.0), 0.001));
  EXPECT_THAT(result.getColor(50, 50), VecNear(Vec3(1.0), 0.001));
}
// Test that the scene looks the same when looket via a mirror
TEST(RenderTest, MirrorScene) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 4.0, kPi / 4.0);
  Scene scene_real(camera);
  scene_real.setEnvironmentLightColor(Vec3(1.0));
  scene_real.addModelFromFile("../models/testing/ball_transparent_diffuse.obj", Vec3(0.0), NormalType::kSmooth);
  scene_real.setSamplingScheme(kImportanceSampling);
  Image image_real = scene_real.render(200, 200, 100, 4, 1);

  Camera camera2(Vec3(0, 0.0, 2.0), Vec3(0.0, 0.0, 0.45),
                Vec3(0.0, 1.0, 0.0), kPi / 4.0, kPi / 4.0);
  Scene scene_mirror(camera2);
  scene_mirror.setEnvironmentLightColor(Vec3(1.0));
  scene_mirror.addModelFromFile("../models/testing/ball_transparent_diffuse.obj", Vec3(0.0), NormalType::kSmooth);
  scene_mirror.addModelFromFile("../models/testing/XY_mirror.obj", Vec3(0.0, 0.0, 2.5), NormalType::kSmooth);
  scene_mirror.setSamplingScheme(kImportanceSampling);
  Image image_mirror = scene_mirror.render(200, 200, 100, 4, 1);
  std::cerr << "Distance between real and mirrored: "
            << image_real.distanceTo(image_mirror) << std::endl;
  EXPECT_LE(image_real.distanceTo(image_mirror), 4.0);
}
// Test that the mirror color is applied correctly
TEST(RenderTest, ColorMirrorEnvironment) {
  Camera camera(Vec3(0, 0.0, 3.0), Vec3(0.0, 0.0, -1.0),
                Vec3(0.0, 1.0, 0.0), kPi / 1.1, kPi / 1.1);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(1.0));
  scene.addModelFromFile("../models/testing/XY_mirror_blue.obj", Vec3(0.0), NormalType::kSmooth);
  scene.setSamplingScheme(kImportanceSampling);
  Image result = scene.render(100, 100, 100, 4, 1);
  EXPECT_THAT(result.getColor(0, 0), VecNear(Vec3(0.0, 0.0, 1.0), 0.001));
  EXPECT_THAT(result.getColor(99, 99), VecNear(Vec3(0.0, 0.0, 1.0), 0.001));
  EXPECT_THAT(result.getColor(20, 36), VecNear(Vec3(0.0, 0.0, 1.0), 0.001));
  EXPECT_THAT(result.getColor(50, 50), VecNear(Vec3(0.0, 0.0, 1.0), 0.001));
}
}  // namespace
