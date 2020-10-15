#include "../../geometry.h"

#include <gmock/gmock-more-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <random>

#include "../../kd_tree.h"
#include "../test_utils.h"

using ::testing::PrintToString;

namespace {
MATCHER_P(VecEq, v, "should equal " + PrintToString(v)) {
  return (v - arg).norm() < EPS;
}
TEST(Vec2Test, ConstructorAllValues) {
  Vec2 v(0.1, -1.1);
  EXPECT_EQ(v[0], 0.1);
  EXPECT_EQ(v[1], -1.1);
}
TEST(Vec2Test, ConstructorOneValue) {
  Vec2 v(1.2345);
  EXPECT_EQ(v[0], 1.2345);
  EXPECT_EQ(v[1], 1.2345);
}
TEST(Vec2Test, DotProduct) {
  Vec2 v1(0.0);
  Vec2 v2(1.1, 2.2);
  EXPECT_NEAR(v1.dot(v2), 0.0, EPS);
  v1 = Vec2(1.1, 2.2);
  v2 = Vec2(0.0, 1.1);
  EXPECT_NEAR(v1.dot(v2), 2.2 * 1.1, EPS);
}
TEST(Vec2Test, SumOfElements) {
  Vec2 v1(0.0);
  EXPECT_NEAR(v1.sum(), 0.0, EPS);
  v1 = Vec2(1.1, -3.4);
  EXPECT_NEAR(v1.sum(), 1.1 - 3.4, EPS);
}
TEST(Vec2Test, MultiplyElementwise) {
  Vec2 v1(0.0);
  Vec2 v2(1.1, 2.2);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec2(0.0, 0.0)));
  v1 = Vec2(1.1, 2.2);
  v2 = Vec2(0.0, 1.1);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec2(0.0, 2.2 * 1.1)));

  v1 = Vec2(1.1, 2.2);
  v2 = Vec2(3.3, 4.4);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec2(1.1 * 3.3, 2.2 * 4.4)));
}
TEST(Vec2Test, PlusAssign) {
  Vec2 v1(1.1, 2.2);
  Vec2 v2(1.1, 1.2);
  v2 += v1;
  EXPECT_NEAR(v2[0], 2.2, EPS);
  EXPECT_NEAR(v2[1], 3.4, EPS);
}
TEST(Vec2Test, Plus) {
  Vec2 v1(1.1, 2.2);
  Vec2 v2(1.1, 1.2);
  Vec2 v3 = v1 + v2;
  EXPECT_NEAR(v3[0], 2.2, EPS);
  EXPECT_NEAR(v3[1], 3.4, EPS);
}
TEST(Vec2Test, MinusAssign) {
  Vec2 v1(1.1, 2.2);
  Vec2 v2(1.1, 1.2);
  v1 -= v2;
  EXPECT_NEAR(v1[0], 0.0, EPS);
  EXPECT_NEAR(v1[1], 1.0, EPS);
}
TEST(Vec2Test, Minus) {
  Vec2 v1(1.1, 2.2);
  Vec2 v2(1.1, 1.2);
  Vec2 v3 = v1 - v2;
  EXPECT_NEAR(v3[0], 0.0, EPS);
  EXPECT_NEAR(v3[1], 1.0, EPS);
}
TEST(Vec2Test, MultipyAssign) {
  Vec2 v1(1.0, 2.2);
  v1 *= 1.1;
  EXPECT_NEAR(v1[0], 1.1, EPS);
  EXPECT_NEAR(v1[1], 2.2 * 1.1, EPS);
}
TEST(Vec2Test, Multipy) {
  Vec2 v1(1.0, 2.2);
  Vec2 v2 = v1 * 1.1;
  EXPECT_NEAR(v2[0], 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 * 1.1, EPS);
  v2 = 1.1 * v1;
  EXPECT_NEAR(v2[0], 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 * 1.1, EPS);
}
TEST(Vec2Test, DivideAssign) {
  Vec2 v1(1.0, 2.2);
  v1 /= 1.1;
  EXPECT_NEAR(v1[0], 1.0 / 1.1, EPS);
  EXPECT_NEAR(v1[1], 2.2 / 1.1, EPS);
}
TEST(Vec2Test, Divide) {
  Vec2 v1(1.0, 2.2);
  Vec2 v2 = v1 / 1.1;
  EXPECT_NEAR(v2[0], 1.0 / 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 / 1.1, EPS);
}
TEST(Vec2Test, Norm) {
  Vec2 v(0.0);
  EXPECT_NEAR(v.norm(), 0.0, EPS);
  v = Vec2(1.0, 0.0);
  EXPECT_NEAR(v.norm(), 1.0, EPS);
  v = Vec2(-1.0, 1.0);
  EXPECT_NEAR(v.norm(), std::sqrt(2), EPS);
  v = Vec2(-13.2, 32.1);
  EXPECT_NEAR(v.norm(), std::sqrt(13.2 * 13.2 + 32.1 * 32.1), EPS);
}
TEST(Vec3Test, ConstructorAllValues) {
  Vec3 v(0.1, -1.1, -2.1);
  EXPECT_EQ(v[0], 0.1);
  EXPECT_EQ(v[1], -1.1);
  EXPECT_EQ(v[2], -2.1);
}
TEST(Vec3Test, ConstructorOneValue) {
  Vec3 v(-1.2345);
  EXPECT_EQ(v[0], -1.2345);
  EXPECT_EQ(v[1], -1.2345);
  EXPECT_EQ(v[2], -1.2345);
}
TEST(Vec3Test, TestDotProduct) {
  Vec3 v1(0.0);
  Vec3 v2(1.1, 2.2, 3.3);
  EXPECT_NEAR(v1.dot(v2), 0.0, EPS);
  v1 = Vec3(1.1, 2.2, 3.3);
  v2 = Vec3(0.0, 1.1, 2.2);
  EXPECT_NEAR(v1.dot(v2), 2.2 * 1.1 + 3.3 * 2.2, EPS);
}
TEST(Vec3Test, SumOfElements) {
  Vec3 v1(0.0);
  EXPECT_NEAR(v1.sum(), 0.0, EPS);
  v1 = Vec3(1.1, 3.4, -4.2);
  EXPECT_NEAR(v1.sum(), 1.1 + 3.4 - 4.2, EPS);
}
TEST(Vec3Test, MultiplyElementwise) {
  Vec3 v1(0.0);
  Vec3 v2(1.1, 2.2, 3.3);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec3(0.0)));
  v1 = Vec3(1.1, 2.2, 3.3);
  v2 = Vec3(0.0, 1.1, 2.2);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec3(0.0, 2.2 * 1.1, 3.3 * 2.2)));
  v1 = Vec3(1.1, 2.2, 3.3);
  v2 = Vec3(4.4, 5.5, 6.6);
  EXPECT_THAT(v1.multiply(v2), VecEq(Vec3(1.1 * 4.4, 2.2 * 5.5, 3.3 * 6.6)));
}
TEST(Vec3Test, PlusAssign) {
  Vec3 v1(1.1, 2.2, 3.2);
  Vec3 v2(1.1, 1.2, -2.1);
  v2 += v1;
  EXPECT_NEAR(v2[0], 2.2, EPS);
  EXPECT_NEAR(v2[1], 3.4, EPS);
  EXPECT_NEAR(v2[2], 1.1, EPS);
}
TEST(Vec3Test, Plus) {
  Vec3 v1(1.1, 2.2, 3.2);
  Vec3 v2(1.1, 1.2, -2.1);
  Vec3 v3 = v1 + v2;
  EXPECT_NEAR(v3[0], 2.2, EPS);
  EXPECT_NEAR(v3[1], 3.4, EPS);
  EXPECT_NEAR(v3[2], 1.1, EPS);
}
TEST(Vec3Test, MinusAssign) {
  Vec3 v1(1.1, 2.2, 3.3);
  Vec3 v2(1.1, 1.2, 3.0);
  v1 -= v2;
  EXPECT_NEAR(v1[0], 0.0, EPS);
  EXPECT_NEAR(v1[1], 1.0, EPS);
  EXPECT_NEAR(v1[2], 0.3, EPS);
}
TEST(Vec3Test, Minus) {
  Vec3 v1(1.1, 2.2, 3.3);
  Vec3 v2(1.1, 1.2, 3.0);
  Vec3 v3 = v1 - v2;
  EXPECT_NEAR(v3[0], 0.0, EPS);
  EXPECT_NEAR(v3[1], 1.0, EPS);
  EXPECT_NEAR(v3[2], 0.3, EPS);
}
TEST(Vec3Test, MultipyAssign) {
  Vec3 v1(1.0, 2.2, 3.3);
  v1 *= 1.1;
  EXPECT_NEAR(v1[0], 1.1, EPS);
  EXPECT_NEAR(v1[1], 2.2 * 1.1, EPS);
  EXPECT_NEAR(v1[2], 3.3 * 1.1, EPS);
}
TEST(Vec3Test, Multipy) {
  Vec3 v1(1.0, 2.2, 3.3);
  Vec3 v2 = v1 * 1.1;
  EXPECT_NEAR(v2[0], 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 * 1.1, EPS);
  EXPECT_NEAR(v2[2], 3.3 * 1.1, EPS);
}
TEST(Vec3Test, DivideAssign) {
  Vec3 v1(1.0, 2.2, 3.3);
  v1 /= 1.1;
  EXPECT_NEAR(v1[0], 1.0 / 1.1, EPS);
  EXPECT_NEAR(v1[1], 2.2 / 1.1, EPS);
  EXPECT_NEAR(v1[2], 3.3 / 1.1, EPS);
}
TEST(Vec3Test, Divide) {
  Vec3 v1(1.0, 2.2, 3.3);
  Vec3 v2 = v1 / 1.1;
  EXPECT_NEAR(v2[0], 1.0 / 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 / 1.1, EPS);
  EXPECT_NEAR(v2[2], 3.3 / 1.1, EPS);
}
TEST(Vec3Test, Norm) {
  Vec3 v(0.0);
  EXPECT_NEAR(v.norm(), 0.0, EPS);
  v = Vec3(1.0, 0.0, 0.0);
  EXPECT_NEAR(v.norm(), 1.0, EPS);
  v = Vec3(-1.0, 1.0, 1.0);
  EXPECT_NEAR(v.norm(), std::sqrt(3), EPS);
  v = Vec3(-13.2, 32.1, 1.0);
  EXPECT_NEAR(v.norm(), std::sqrt(13.2 * 13.2 + 32.1 * 32.1 + 1.0), EPS);
}
TEST(OnSameSideOfPlane, Simple) {
  Vec3 a(1.0, 1.0, 0.0);
  Vec3 b(0.0, 1.0, 1.0);
  Vec3 n(0.0, 1.0, 0.0);
  EXPECT_TRUE(onSameSideOfPlane(a, b, n));
  a = Vec3(1.0, 1.0, 0.0);
  b = Vec3(0.0, -1.0, 1.0);
  n = Vec3(0.0, 1.0, 0.0);
  EXPECT_FALSE(onSameSideOfPlane(a, b, n));
}
TEST(OnSameSideOfPlane, IdenticalVectors) {
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(1.0, 2.0, 3.0);
  Vec3 n(0.0, 1.0, 0.0);
  EXPECT_TRUE(onSameSideOfPlane(a, b, n));
}
TEST(Project, IdenticalVectors) {
  Vec3 a(1.2, 2.3, 3.4);
  Vec3 b(1.2, 2.3, 3.4);
  EXPECT_THAT(project(a, b), VecEq(Vec3(1.2, 2.3, 3.4)));
}
TEST(Project, Simple) {
  Vec3 a(1.2, 2.3, 3.4);
  Vec3 b(-12.0, 0.0, 0.0);
  EXPECT_THAT(project(a, b), VecEq(Vec3(1.2, 0.0, 0.0)));
  a = Vec3(1.2, 2.3, 3.4);
  b = Vec3(-12.0, 0.0, 0.0);
  EXPECT_THAT(project(a, b), VecEq(Vec3(1.2, 0.0, 0.0)));
}
TEST(MirrorOver, IdenticaVectors) {
  Vec3 a(1.2, 2.3, 3.4);
  Vec3 b(1.2, 2.3, 3.4);
  EXPECT_THAT(mirrorOver(a, b), VecEq(Vec3(1.2, 2.3, 3.4)));
}
TEST(MirrorOver, Simple) {
  Vec3 a(1.1, 0.0, 0.0);
  Vec3 b(1.0, 1.0, 0.0);
  EXPECT_THAT(mirrorOver(a, b), VecEq(Vec3(0.0, 1.1, 0.0)));
}
TEST(MirrorOver, OverNormal) {
  Vec3 a(1.1, 2.1, -2.3);
  Vec3 b(-2.1, 1.1, 0.0);
  EXPECT_THAT(mirrorOver(a, b), VecEq(-1.0 * a));
}
// Checks that sum of a random vector and the mirror image lie on
// correct vector and that both are of correct length
TEST(MirrorOver, Random) {
  std::mt19937 mt(1337);
  for(int i = 0; i < 10000; ++i) {
    Vec3 a = test::randomVec3(0.0, 10.0, mt);
    Vec3 b = test::randomVec3(0.0, 10.0, mt);
    Vec3 mirrored = mirrorOver(a, b);
    ASSERT_NEAR(a.norm(), mirrored.norm(), EPS);
    Vec3 sum = a + mirrored;
    sum /= sum.norm();
    b /= b.norm();
    ASSERT_NEAR(sum.dot(b), 1.0, EPS);
  }
  Vec3 a(1.1, 2.1, -2.3);
  Vec3 b(-2.1, 1.1, 0.0);
  EXPECT_THAT(mirrorOver(a, b), VecEq(-1.0 * a));
}
TEST(VoxelTest, IntersectsSimple) {
  Voxel v({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  Ray r({2.0, 2.0, 2.0}, {1.0, 1.0, 1.0});
  EXPECT_FALSE(v.intersects(r));
  r = Ray({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0});
  EXPECT_TRUE(v.intersects(r));
  r = Ray({2.0, 2.0, 2.0}, {-1.0, -1.0, -1.0});
  EXPECT_TRUE(v.intersects(r));
}
// Random test for rays intersecting inner points of the voxel
TEST(VoxelTest, IntersectsRandom) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    std::uniform_real_distribution<double> dist(-1, 1);
    Voxel v(test::randomVec3(-1, 1, mt), test::randomVec3(-1, 1, mt));

    Vec3 inner_point;
    for (int j = 0; j < 3; ++j) {
      if (v.lo[j] > v.hi[j]) std::swap(v.lo[j], v.hi[j]);
      std::uniform_real_distribution<double> inner_dist(v.lo[j], v.hi[j]);
      inner_point[j] = inner_dist(mt);
    }
    bool ok = 1;
    for (int j = 0; j < 3; ++j) {
      // make sure inner point doesn't lie on surface
      if (inner_point[j] < v.lo[j] + EPS || inner_point[j] > v.hi[j] - EPS) {
        ok = 0;
      }
    }
    if (ok) {
      Vec3 other_point = test::randomVec3(-1, 1, mt);
      double scale = test::randomLogUniformReal(-10, 10, mt);
      Vec3 direction = scale * (inner_point - other_point);
      Ray r(other_point, direction);
      ASSERT_TRUE(v.intersects(r)) << "v = " << v << "\nr = " << r;
    }
  }
}
// TODO better checks for non intersecting rays
// TODO better check for axis-aligned rays
TEST(VoxelTest, Clip) {
  Voxel v({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  AxisPlane plane({0, 0.5});
  v.clip(plane, 1);
  Voxel c({0.5, 0.0, 0.0}, {1.0, 1.0, 1.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);

  v = Voxel({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  plane = AxisPlane({2, 0.5});
  v.clip(plane, 0);
  c = Voxel({0.0}, {1.0, 1.0, 0.5});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);
}
TEST(VoxelTest, CoverPoint) {
  Voxel v({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  Vec3 p(0.0, 0.0, 0.0);
  v.cover(p);
  Voxel c({0.0}, {0.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);

  v = Voxel({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  p = Vec3(0.1, 1.2, 2.3);
  v.cover(p);
  c = Voxel({0.0}, p);
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);

  v = Voxel({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  p = Vec3(-0.1, 1.2, -2.3);
  v.cover(p);
  c = Voxel({-0.1, 0.0, -2.3}, {0.0, 1.2, 0.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);

  v = Voxel({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  p = Vec3(0.5, 1.0, -0.0);
  v.cover(p);
  c = Voxel({0.0}, {1.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);
}
TEST(VoxelTest, CoverTriangle) {
  Voxel v({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  Triangle t{{0.0, -1.0, 0.0}, {-2.2, 0.0, 1.0}, {5.2, 0.0, -1.2}};
  v.cover(t);
  Voxel c({-2.2, -1.0, -1.2}, {5.2, 0.0, 1.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);
}
TEST(VoxelTest, Area) {
  Voxel v({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  EXPECT_NEAR(v.area(), 0, EPS);
  v = Voxel({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
  EXPECT_NEAR(v.area(), 0, EPS);
  v = Voxel({0.0, 0.0, 0.0}, {1.0, 1.0, 0.0});
  EXPECT_NEAR(v.area(), 2.0, EPS);
  v = Voxel({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  EXPECT_NEAR(v.area(), 6.0, EPS);
}
TEST(TriangleTest, PointFromBarySimple) {
  Triangle t({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  Vec2 bary_coords(0.0, 0.0);
  Vec3 r(t.pointFromBary(bary_coords));
  Vec3 c(0.0, 0.0, 0.0);
  EXPECT_NEAR((r - c).norm(), 0, EPS);

  t = Triangle({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  bary_coords = Vec2(0.2, 0.43);
  r = t.pointFromBary(bary_coords);
  c = Vec3(0.2, 0.43, 0.0);
  EXPECT_NEAR((r - c).norm(), 0, EPS);
}
// tests that random barycentric coordinates will be on the triangle plane
TEST(TriangleTest, PointFromBaryRandom) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    Triangle t(test::randomVec3(-1.0, 1.0, mt), test::randomVec3(-1.0, 1.0, mt),
               test::randomVec3(-1.0, 1.0, mt));
    Vec2 bary_coords(test::randomBaryCoords(mt));
    Vec3 p(t.pointFromBary(bary_coords));
    Vec3 dot_normal = (p - t.p0).cross(t.p1 - t.p0).dot(t.p2 - t.p0);
    ASSERT_NEAR(dot_normal.norm(), 0, EPS);
  }
}
TEST(TriangleTest, Area) {
  Triangle t({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  EXPECT_NEAR(t.area(), 0.0, EPS);
  t = Triangle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
  EXPECT_NEAR(t.area(), 0.0, EPS);
  t = Triangle({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  EXPECT_NEAR(t.area(), 0.5, EPS);
  t = Triangle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 1.0});
  EXPECT_NEAR(t.area(), std::sqrt(1.0 / 2.0), EPS);
}
TEST(TriangleTest, Translate) {
  Triangle t(Vec3(0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  t.translate(Vec3(1.0, 2.0, 3.0));
  EXPECT_THAT(t.p0, VecEq(Vec3(1.0, 2.0, 3.0)));
  EXPECT_THAT(t.p1, VecEq(Vec3(2.0, 2.0, 3.0)));
  EXPECT_THAT(t.p2, VecEq(Vec3(1.0, 3.0, 3.0)));
}
TEST(TriangleTest, RayIntersectionSimple) {
  Triangle t({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  Ray r({0.2, 0.2, 0.2}, {0.0, 0.0, 1.0});
  RayIntersection ri = t.getRayIntersection(r);
  EXPECT_NEAR(ri.distance, -0.2, EPS);

  t = Triangle({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  r = Ray({2.0, 2.0, 2.0}, {0.0, 0.0, 1.0});
  ri = t.getRayIntersection(r);
  EXPECT_EQ(ri.distance, INF);
}
TEST(TriangleTest, RayIntersectionParallel) {
  Triangle t({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  Ray r({0.0, 0.0, 0.0}, {1.0, 1.0, 0.0});
  RayIntersection ri = t.getRayIntersection(r);
  EXPECT_EQ(ri.distance, INF);
}
// test that getRayIntersection identifies all intersections correctly
TEST(TriangleTest, RayIntersectionRandomIntersects) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    Triangle random_t(test::randomVec3(-1.0, 1.0, mt),
                      test::randomVec3(-1.0, 1.0, mt),
                      test::randomVec3(-1.0, 1.0, mt));

    Vec2 bary_coords = test::randomBaryCoords(mt);
    // skip coordinates that are too close to border
    if (bary_coords[0] < EPS || bary_coords[1] < EPS ||
        bary_coords[0] + bary_coords[1] > 1 - EPS) {
      continue;
    }
    Vec3 inner_point(random_t.pointFromBary(bary_coords));
    Vec3 ray_origin(test::randomVec3(-1.0, 1.0, mt));
    // ray origin should not be on the plane
    // TODO implement this to geometry.cpp?
    Vec3 tmp = (ray_origin - random_t.p0).cross(random_t.p1 - random_t.p0);
    if (std::abs(tmp.dot(random_t.p2 - random_t.p0)) < EPS) {
      continue;
    }
    double scale = test::randomLogUniformReal(-1.0, 1.0, mt);
    // also test negative scale (ray going away from the triangle)
    if (i % 2 == 1) {
      scale *= -1;
    }
    Ray r(ray_origin, scale * (inner_point - ray_origin));
    RayIntersection ri = random_t.getRayIntersection(r);
    ASSERT_NEAR(ri.distance, 1 / scale, EPS)
        << "r = " << r << "\nt = " << random_t << "\nbary = " << bary_coords
        << std::endl;
  }
}
// test that getRayIntersection identifies all non-intersections correctly
TEST(TriangleTest, RayIntersectionRandomDoesNotIntersect) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    Triangle random_t(test::randomVec3(-1.0, 1.0, mt),
                      test::randomVec3(-1.0, 1.0, mt),
                      test::randomVec3(-1.0, 1.0, mt));

    // generate barycentric coordinates that can be outside the triangle
    std::uniform_real_distribution bary_distribution(-100.0, 100.0);
    Vec2 bary_coords(bary_distribution(mt), bary_distribution(mt));
    // skip coordinates thaat are inside the triangle
    if (bary_coords[0] > -TEST_EPS && bary_coords[0] < 1 + TEST_EPS &&
        bary_coords[1] > -TEST_EPS && bary_coords[1] < 1 + TEST_EPS &&
        bary_coords[0] + bary_coords[1] > -TEST_EPS &&
        bary_coords[0] + bary_coords[1] < 1 + TEST_EPS) {
      continue;
    }
    // pointFromBary only supports points on the triangle
    Vec3 out_point = random_t.p0 * (1 - bary_coords[0] - bary_coords[1]) +
                     random_t.p1 * bary_coords[0] +
                     random_t.p2 * bary_coords[1];
    Vec3 ray_origin(test::randomVec3(-1.0, 1.0, mt));
    // ray origin should not be on the plane
    Vec3 tmp = (ray_origin - random_t.p0).cross(random_t.p1 - random_t.p0);
    if (std::abs(tmp.dot(random_t.p2 - random_t.p0)) < EPS) {
      continue;
    }
    double scale = test::randomLogUniformReal(-1.0, 1.0, mt);
    // also test negative scale (ray going away from the triangle)
    if (i % 2 == 1) {
      scale *= -1;
    }
    Ray r(ray_origin, scale * (out_point - ray_origin));
    RayIntersection ri = random_t.getRayIntersection(r);
    ASSERT_EQ(ri.distance, INF) << "r = " << r << "\nt = " << random_t
                                << "\nbary = " << bary_coords << std::endl;
  }
}
TEST(PlanePolygon, ConstructorPoints) {
  PlanePolygon p(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  std::vector<Vec3> r = {p[0], p[1], p[2]};
  EXPECT_NEAR((p[0] - Vec3(1.0, 0.0, 0.0)).norm(), 0, EPS);
  EXPECT_NEAR((p[1] - Vec3(1.1)).norm(), 0, EPS);
  EXPECT_NEAR((p[2] - Vec3(2.2)).norm(), 0, EPS);
}
TEST(PlanePolygon, ConstructorTriangle) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(t);
  std::vector<Vec3> r = {p[0], p[1], p[2]};
  EXPECT_NEAR((p[0] - Vec3(1.0, 0.0, 0.0)).norm(), 0, EPS);
  EXPECT_NEAR((p[1] - Vec3(1.1)).norm(), 0, EPS);
  EXPECT_NEAR((p[2] - Vec3(2.2)).norm(), 0, EPS);
}
TEST(PlanePolygon, BoundingBoxSimple) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(t);
  Voxel r(p.getBoundingBox());
  Voxel c(Vec3(1.0, 0.0, 0.0), Vec3(2.2));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  p = PlanePolygon(t);
  r = Voxel(p.getBoundingBox());
  c = Voxel(Vec3(-12.23, -15.2, 2.2), Vec3(1.1, 2.2, 3.9999));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);
}
TEST(PlanePolygon, IntersectSimple) {
  Triangle t(Vec3(0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 1.0));
  PlanePolygon p(t);
  p.intersect({2, 0.5}, 0);
  Voxel r = p.getBoundingBox();
  Voxel c = Voxel(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.5, 0.5));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);
}
// TODO add more intersection tests
TEST(PlanePolygon, Size) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(t);
  EXPECT_EQ(p.size(), 3);
}
TEST(ClipTriangle, MinMin) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(t);
  EXPECT_NEAR(ct.min(0), -12.23, EPS);
  EXPECT_NEAR(ct.min(1), -15.2, EPS);
  EXPECT_NEAR(ct.min(2), 2.2, EPS);
  EXPECT_NEAR(ct.max(0), 1.1, EPS);
  EXPECT_NEAR(ct.max(1), 2.2, EPS);
  EXPECT_NEAR(ct.max(2), 3.9999, EPS);
}
TEST(ClipTriangle, IsAxisAligned) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  t = Triangle(Vec3(1.0, 0.0, 2.5), Vec3(1.0, 2.2, 3.9999),
               Vec3(1.0, -15.2, 2.2));
  ct = ClipTriangle(t);
  EXPECT_TRUE(ct.isAxisAligned(0));
  EXPECT_FALSE(ct.isAxisAligned(1));
  EXPECT_FALSE(ct.isAxisAligned(2));

  t = Triangle(Vec3(1.0, 2.0, 12), Vec3(2.0, 2.0, 3.9999), Vec3(3.0, 2.0, 2.2));
  ct = ClipTriangle(t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  EXPECT_TRUE(ct.isAxisAligned(1));
  EXPECT_FALSE(ct.isAxisAligned(2));

  t = Triangle(Vec3(1.0, 12, 3.3), Vec3(2.0, 2.0, 3.3), Vec3(3.0, 3.0, 3.3));
  ct = ClipTriangle(t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  EXPECT_FALSE(ct.isAxisAligned(1));
  EXPECT_TRUE(ct.isAxisAligned(2));
}
TEST(ClipTriangle, OverlapSidesSimple) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(t);
  std::pair<bool, bool> r = ct.overlapsSides({1, 1.0}, 1);
  EXPECT_EQ(r, std::make_pair(true, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({0, -100.0}, 1);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({2, 100.0}, 0);
  EXPECT_EQ(r, std::make_pair(true, false));
}
TEST(ClipTriangle, OverlapSidesParallel) {
  Triangle t(Vec3(0.0, 0.0, 2.5), Vec3(0.0, 2.2, 3.9999),
             Vec3(0.0, -15.2, 2.2));
  ClipTriangle ct(t);
  std::pair<bool, bool> r = ct.overlapsSides({0, 0.0}, 1);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(1.0, 1.1, 2.5), Vec3(1.2, 1.1, 3.9999),
               Vec3(1.1, 1.1, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({1, 1.1}, 0);
  EXPECT_EQ(r, std::make_pair(true, false));
}
// tests the case where the plane intersects a zero area part of the
// non axis-aligned triangle
TEST(ClipTriangle, OverlapSidesZeroArea) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(t);
  std::pair<bool, bool> r = ct.overlapsSides({0, -12.23}, 1);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({0, -12.23}, 0);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({1, 2.2}, 0);
  EXPECT_EQ(r, std::make_pair(true, false));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(t);
  r = ct.overlapsSides({1, 2.2}, 1);
  EXPECT_EQ(r, std::make_pair(true, false));
}
TEST(ClipTriangle, ClipSimple) {
  Triangle t(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ClipTriangle ct(t);
  ct.clip({0, 0.5}, 0);
  EXPECT_NEAR(ct.max(0), 0.5, EPS);
  EXPECT_NEAR(ct.max(1), 1.0, EPS);
  EXPECT_NEAR(ct.min(0), 0.0, EPS);

  t = Triangle(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ct = ClipTriangle(t);
  ct.clip({0, 0.5}, 1);
  EXPECT_NEAR(ct.min(0), 0.5, EPS);
  EXPECT_NEAR(ct.max(1), 0.5, EPS);
}
TEST(ClipTriangle, ClipParallel) {
  Triangle t(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ClipTriangle ct(t);
  ct.clip({2, 0.0}, 0);
  EXPECT_NEAR(ct.min(2), 0.0, EPS);
  EXPECT_NEAR(ct.max(2), 0.0, EPS);
  EXPECT_NEAR(ct.min(0), 0.0, EPS);
  EXPECT_NEAR(ct.max(0), 1.0, EPS);
}
TEST(FirstRayTriangleIntersection, Simple) {
  std::vector<Triangle> scene;
  scene.push_back(
      Triangle(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)));
  scene.push_back(
      Triangle(Vec3(2.0, 2.0, 0.0), Vec3(2.5, 2.0, 0.0), Vec3(2.0, 2.5, 0.0)));
  scene.push_back(
      Triangle(Vec3(2.2, 1.0, 0.0), Vec3(1.3, 2.2, 0.0), Vec3(1.1, 2.1, 0.0)));
  Ray r = Ray(Vec3(0.0, 0.0, 1.0), Vec3(2.1, 2.1, 1.0));
  RayTriangleIntersection rti = firstRayTriangleIntersection(scene, r);
  EXPECT_EQ(rti.index, 3);

  r = Ray(Vec3(0.0, 0.0, 1.0), Vec3(2.1, 2.1, -1.0));
  rti = firstRayTriangleIntersection(scene, r);
  EXPECT_EQ(rti.index, 1);

  scene.clear();
  scene.push_back(
      Triangle(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)));
  scene.push_back(
      Triangle(Vec3(0.0, 0.0, 1.0), Vec3(1.0, 0.0, 1.0), Vec3(0.0, 1.0, 1.0)));
  scene.push_back(
      Triangle(Vec3(0.0, 0.0, 2.0), Vec3(1.0, 0.0, 2.0), Vec3(0.0, 1.0, 2.0)));

  r = Ray(Vec3(0.3, 0.3, -10.0), Vec3(0.0, 0.0, 1.0));
  rti = firstRayTriangleIntersection(scene, r);
  EXPECT_EQ(rti.index, 0);
  EXPECT_THAT(rti.bary_coords, VecEq(Vec2(0.3, 0.3)));

  r = Ray(Vec3(0.2, 0.2, 0.5), Vec3(0.2, 0.2, 1.0));
  rti = firstRayTriangleIntersection(scene, r);
  EXPECT_EQ(rti.index, 1);
  EXPECT_THAT(rti.bary_coords, VecEq(Vec2(0.3, 0.3)));

  r = Ray(Vec3(0.2, 0.2, 2 - 0.01), Vec3(0.0, 0.0, 1.0));
  rti = firstRayTriangleIntersection(scene, r);
  EXPECT_EQ(rti.index, 2);
  EXPECT_THAT(rti.bary_coords, VecEq(Vec2(0.2, 0.2)));
}
// Tests that the function always finds an intersection that
// is at least as close as the correct intersection
// Also tests that the intersections are really on the ray
TEST(FirstRayTriangleIntersection, Random) {
  std::mt19937 mt(1337);
  int n_tests_run = 0;
  int n_same_triangle = 0;
  for (int i = 0; i < 10'000; ++i) {
    double max_triangle_size = test::randomLogUniformReal(-10, 0, mt);
    std::vector<Triangle> scene =
        test::randomTriangleVector(-1, 1, max_triangle_size, 100, mt);
    Vec2 bary_coords = test::randomBaryCoords(mt);
    Vec3 p = scene[0].pointFromBary(bary_coords);
    Vec3 ray_origin = p + test::randomVec3(-0.04, 0.04, mt);
    // 0.001 should always still be a lot larger than EPS
    if ((ray_origin - p).norm() < 0.001) continue;
    if (test::pointOnTrianglePlane(scene[0], ray_origin)) continue;

    double scale = test::randomLogUniformReal(-4, 10, mt);
    Ray r(ray_origin, scale * (p - ray_origin));
    RayTriangleIntersection rti = firstRayTriangleIntersection(scene, r);
    ASSERT_NE(rti.index, scene.size());
    Vec3 intersection = scene[rti.index].pointFromBary(rti.bary_coords);
    // check that the intersection found is not
    // further away than it should be
    ASSERT_LE((intersection - ray_origin).norm(),
              (p - ray_origin).norm() + EPS);
    // check that the intersection is on the ray
    ASSERT_NEAR((intersection - ray_origin).cross(p - ray_origin).norm(), 0,
                EPS);
    if (rti.index == 0) ++n_same_triangle;
    ++n_tests_run;
  }
  ASSERT_GE(n_tests_run, 7000) << "Problem in generating the test cases";
  std::cerr << "random tests run: " << n_tests_run << std::endl;
  std::cerr << "of these, " << n_same_triangle
            << " intersected first the triangle used to generate the ray "
            << std::endl;
}
TEST(PointOnSegment, Simple) {
  Vec3 a(0.0);
  Vec3 b(1.0);
  Vec3 c(2.0);
  EXPECT_TRUE(pointOnSegment(b, a, c));
  EXPECT_FALSE(pointOnSegment(a, b, c));
  a = Vec3(0.0);
  b = Vec3(0.0, 1.0, 1.0);
  c = Vec3(0.0, 0.0, 2.0);
  EXPECT_FALSE(pointOnSegment(b, a, c));
}
// This also tests the current implementation special case b[0] == 0
TEST(RotateYTo, Simple) {
  Vec3 a(0.0, 1.0, 0.0);
  Vec3 b(0.0, 0.0, 1.0);
  EXPECT_THAT(rotateYTo(a, b), VecEq(b));
  b = Vec3(1.0, 0.0, 0.0);
  EXPECT_THAT(rotateYTo(a, b), VecEq(b));
  b = Vec3(0.0, 1.0, 0.0);
  EXPECT_THAT(rotateYTo(a, b), VecEq(b));
  b = Vec3(0.0, 0.0, 1.0);
  EXPECT_THAT(rotateYTo(a, b), VecEq(b));
  b = Vec3(1.0, 2.0, 1.23);
  b /= b.norm();
  EXPECT_THAT(rotateYTo(a, b), VecEq(b));
}
TEST(UniformRandomHemispherePoint, PointsOnHemisphere) {
  std::mt19937 mt(1337);
  Vec3 direction = Vec3(1.1, 2.2, 3.3);
  for(int i = 0; i < 10000; ++i) {
    Vec3 random_point = uniformRandomHemispherePoint(direction, mt);
    ASSERT_GE(random_point.dot(direction), -EPS);
    ASSERT_NEAR(random_point.norm(), 1, EPS);
  }
}
TEST(UniformRandomHemispherePoint, UniformDistribution) {
  std::mt19937 mt(1337);
  Vec3 direction = Vec3(-3.3, -2.2, -1.1);
  Vec3 normalized_direction = direction / direction.norm();
  // direction.dot(direction_normal) == 0
  // if the distribution is uniform, there should be exactly
  // half as many dots wihin some angle from direction_normal
  // than wihtin some angle from direction
  Vec3 direction_normal(direction[1],
                        -direction[0], 0.0);
  direction_normal = direction_normal / direction_normal.norm();
  assert(direction.dot(direction_normal) < EPS);
  Vec3 in_between = direction_normal +  direction;
  in_between = in_between / in_between.norm();
  int n_points = 500000;
  int n_near_direction = 0;
  int n_near_in_between = 0;
  int n_near_normal = 0;
  for (int i = 0; i < n_points; ++i) {
    Vec3 random_point = uniformRandomHemispherePoint(direction, mt);
    // project to the normalized direction
    if (random_point.dot(normalized_direction) > 0.8) {
      ++n_near_direction;
    }
    if (random_point.dot(direction_normal) > 0.8) {
      ++n_near_normal;
    }
    if(random_point.dot(in_between) > 0.8) {
      ++n_near_in_between;
    }
  }
  double fraction_near_direction = (double)n_near_direction / n_points;
  double fraction_near_normal = (double)n_near_normal / n_points;
  double fraction_near_in_between = (double)n_near_in_between / n_points;
  EXPECT_NEAR(fraction_near_normal*2, fraction_near_direction, 0.01);
  EXPECT_NEAR(fraction_near_in_between, fraction_near_direction, 0.01);
}
TEST(CosineExponentRandomPoint, PointsOnHemisphere) {
  std::mt19937 mt(1337);
  Vec3 direction = Vec3(1.1, 2.2, 3.3);
  std::uniform_real_distribution exponent_distribution(0.001, 5.0);
  for(int i = 0; i < 10000; ++i) {
    double exponent = exponent_distribution(mt);
    Vec3 random_point = cosineExponentRandomPoint(direction, exponent, mt);
    ASSERT_GE(random_point.dot(direction), -EPS);
    ASSERT_NEAR(random_point.norm(), 1, EPS);
  }
}
}  // namespace
