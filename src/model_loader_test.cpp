#include "model_loader.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>

#include "geometry.h"
#include "raycaster.h"

using ::testing::ElementsAreArray;
using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
MATCHER_P(TriangleEq, v, "should equal " + PrintToString(v)) {
  if ((v.p0 - arg.p0).norm() > EPS) return 0;
  if ((v.p1 - arg.p1).norm() > EPS) return 0;
  if ((v.p2 - arg.p2).norm() > EPS) return 0;
  return 1;
}
TEST(TestLoadModel, SimpleTriangle) {
  Model model;
  loadModel("../models/green_triangle.obj", model, kRough);
  Triangle c(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f),
             Vec3(0.0f, 1.0f, 0.0f));
  ASSERT_EQ(model.scene_triangles.size(), 1);
  EXPECT_THAT(model.scene_triangles[0].triangle, TriangleEq(c));

  EXPECT_THAT(model.scene_triangles[0].triangle, TriangleEq(c));
}
TEST(TestLoadModel, MaterialDiffuse) {
  Model model;
  loadModel("../models/green_triangle.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material material = model.scene_triangles[0].material;
  EXPECT_THAT(material.diffuse, VecEq(Vec3(0.0f, 0.8f, 0.0f)));
  EXPECT_THAT(material.emitted, VecEq(Vec3(0.0f, 0.0f, 0.0f)));
}
TEST(TestLoadModel, MaterialEmitted) {
  Model model;
  loadModel("../models/white_light.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material material = model.scene_triangles[0].material;
  EXPECT_THAT(material.diffuse, VecEq(Vec3(0.0f, 0.0f, 0.0f)));
  EXPECT_THAT(material.emitted, VecEq(Vec3(1.0f, 1.0f, 1.0f)));
}
TEST(TestLoadModel, MaterialSpecular) {
  Model model;
  loadModel("../models/green_triangle.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material material = model.scene_triangles[0].material;
  EXPECT_THAT(material.specular, VecEq(Vec3(0.2f, 0.2f, 0.2f)));
  EXPECT_NEAR(material.specular_exp, 2.0f, EPS);
}
TEST(TestLoadModel, MaterialTransparent) {
  Model model;
  loadModel("../models/greenish_transparent_triangle.obj", model, kRough);
  loadModel("../models/greenish_transparent_triangle.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material material = model.scene_triangles[0].material;
  EXPECT_THAT(material.transparent, VecEq(Vec3(1.0f, 0.2f, 0.95f)));
}
TEST(TestLoadModel, MaterialIndexOfRefraction) {
  Model model;
  loadModel("../models/greenish_transparent_triangle.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material material = model.scene_triangles[0].material;
  EXPECT_EQ(material.index_of_refraction, 1.5f);
}
TEST(TestLoadModel, Normals) {
  Model model;
  loadModel("../models/green_triangle.obj", model, kRough);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  EXPECT_THAT(model.scene_triangles[0].normals[0],
              VecEq(Vec3(-1.0f, 0.0f, 0.0f)));
  EXPECT_THAT(model.scene_triangles[0].normals[1],
              VecEq(Vec3(-1.0f, 0.0f, 0.0f)));
  EXPECT_THAT(model.scene_triangles[0].normals[2],
              VecEq(Vec3(-1.0f, 0.0f, 0.0f)));
}
}  // namespace
