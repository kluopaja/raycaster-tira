#include "model_loader.h"
#include "raycaster.h"
#include "geometry.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>

using ::testing::PrintToString;
using ::testing::ElementsAreArray;

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
  loadModel("../models/green_triangle.obj", model);
  Triangle c(Vec3(0.0, 0.0, 0.0),
             Vec3(0.0, 0.0, 1.0),
             Vec3(0.0, 1.0, 0.0));
  ASSERT_EQ(model.scene_triangles.size(), 1);
  EXPECT_THAT(model.scene_triangles[0]->triangle, TriangleEq(c));

  EXPECT_THAT(model.scene_triangles[0]->triangle, TriangleEq(c));
}
TEST(TestLoadModel, MaterialDiffuse) {
  Model model;
  loadModel("../models/green_triangle.obj", model);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material* material = model.scene_triangles[0]->material;
  ASSERT_NE(material, nullptr);
  EXPECT_THAT(material->diffuse, VecEq(Vec3(0.0, 1.0, 0.0)));
  EXPECT_THAT(material->emitted, VecEq(Vec3(0.0, 0.0, 0.0)));
}
TEST(TestLoadModel, MaterialEmitted) {
  Model model;
  loadModel("../models/white_light.obj", model);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  Material* material = model.scene_triangles[0]->material;
  ASSERT_NE(material, nullptr);
  EXPECT_THAT(material->diffuse, VecEq(Vec3(0.0, 0.0, 0.0)));
  EXPECT_THAT(material->emitted, VecEq(Vec3(1.0, 1.0, 1.0)));
}
TEST(TestLoadModel, Normals) {
  Model model;
  loadModel("../models/green_triangle.obj", model);
  ASSERT_EQ(model.scene_triangles.size(), 1);
  EXPECT_THAT(model.scene_triangles[0]->normals[0], VecEq(Vec3(-1.0, 0.0, 0.0)));
  EXPECT_THAT(model.scene_triangles[0]->normals[1], VecEq(Vec3(-1.0, 0.0, 0.0)));
  EXPECT_THAT(model.scene_triangles[0]->normals[2], VecEq(Vec3(-1.0, 0.0, 0.0)));
}
} // namespace
