#include "../../input.h"

#include <cmath>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../../geometry.h"
#include "../../vector.h"

using ::testing::PrintToString;
namespace {
bool VecEq(const Vec3& a, const Vec3& b) {
  return (a - b).norm() < EPS;
}
bool DoubleEq(double a, double b) {
  return std::abs(a - b) < EPS;
}
bool operator==(const CameraConfig& a, const CameraConfig& b) {
  return VecEq(a.pos, b.pos) &&
         VecEq(a.front, b.front) &&
         VecEq(a.up, b.up) &&
         DoubleEq(a.x_fov, b.x_fov) &&
         DoubleEq(a.y_fov, b.y_fov);
}
bool operator==(const ModelConfig& a, const ModelConfig& b) {
  return (a.file == b.file) &&
         VecEq(a.pos, b.pos) &&
         a.normal == b.normal;
}
bool operator==(const PointLightConfig& a, const PointLightConfig& b) {
  return VecEq(a.pos, b.pos) &&
         VecEq(a.color, b.color);
}
bool operator==(const EnvironmentLightConfig& a,
                const EnvironmentLightConfig& b) {
  return VecEq(a.color, b.color) && 
         (a.type == b.type) &&
         VecEq(a.direction, b.direction) &&
         DoubleEq(a.exp, b.exp);
}
bool operator==(const RenderConfig& a, const RenderConfig& b) {
  return (a.width == b.width) &&
         (a.height == b.height) &&
         (a.pixel_rays == b.pixel_rays) &&
         (a.depth == b.depth) &&
         (a.branching == b.branching);
}
bool operator==(const ImageConfig& a, const ImageConfig& b) {
  return (a.file == b.file) &&
         DoubleEq(a.truncate, b.truncate) &&
         DoubleEq(a.scale_max, b.scale_max);
}
MATCHER_P(RecipeEq, r, "should equal " + PrintToString(r)) {
  if(!(r.camera == arg.camera)) return 0;
  if(!(r.models.size() == arg.models.size())) return 0;
  if(!(r.models == arg.models)) return 0;
  if(!(r.point_lights == arg.point_lights)) return 0;
  if(!(r.environment_light == arg.environment_light)) return 0;
  if(!(r.render == arg.render)) return 0; 
  if(!(r.images == arg.images)) return 0;
  return 1;
}
TEST(LoadRecipe, Empty) {
  Recipe r = loadRecipe("../recipes/testing/recipe_empty.txt", 0);
  Recipe correct;
  EXPECT_THAT(r, RecipeEq(correct));
}
// Test loading empty sections (only containing the required 'file'
// field in 'model' and 'image'
TEST(LoadRecipe, EmptySections) {
  Recipe r = loadRecipe("../recipes/testing/recipe_empty_sections.txt", 0);
  Recipe correct;
  ModelConfig model;
  model.file = "../recipes/testing/file.txt";
  correct.models.pushBack(model);
  correct.models.pushBack(model);
  correct.point_lights.pushBack({});
  correct.point_lights.pushBack({});
  ImageConfig image;
  image.file = "../recipes/testing/file.txt";
  correct.images.pushBack(image);
  correct.images.pushBack(image);
  EXPECT_THAT(r, RecipeEq(correct));
}
TEST(LoadRecipe, AllDefined) {
  Recipe r = loadRecipe("../recipes/testing/recipe_all_defined_1.txt", 0);
  Recipe correct;
  correct.camera.pos = Vec3(0.0, 0.0, 10.0);
  correct.camera.front = Vec3(0.0, 1.0, -1.0);
  correct.camera.up = Vec3(0.0, 1.0, 1.0);
  correct.camera.x_fov = 90 * kPi / 180.0;
  correct.camera.y_fov = 85 * kPi / 180.0;

  correct.models.pushBack({"../recipes/testing/ball.obj",
                           {10.0, 10.0, 123.0}, NormalType::kRough});
  correct.models.pushBack({"../recipes/testing/triangle.obj",
                           {-1.0, 2.0, 3.0}, NormalType::kSmooth});
  correct.point_lights.pushBack({{1.0, 2.0, 3.0}, {100.0, 100.0, 100.0}});
  correct.point_lights.pushBack({{3.0, 2.0, 1.0}, {4.0, 5.0, 6.0}});

  correct.environment_light.color = Vec3(1.0, 1.0, 1.0);
  correct.environment_light.type = EnvironmentLightType::kDirected;
  correct.environment_light.direction = Vec3(10.0, 10.0, 10.0);
  correct.environment_light.exp = 1000.0;

  correct.render.width = 100;
  correct.render.height = 200;
  correct.render.pixel_rays = 123;
  correct.render.depth = 10;
  correct.render.branching = 14;

  correct.images.pushBack({"asdf.ppm", 0.99, 0.9});
  correct.images.pushBack({"asdf.ppm", 0.89, 0.19});
  correct.render.depth = 10;
}

}   // namespace
