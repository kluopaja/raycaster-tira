#include "geometry.hpp"

#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>

#include <random>
using ::testing::UnorderedElementsAre;
namespace {

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
TEST(Vec2Test, Minus) {
  Vec2 v1(1.1, 2.2);
  Vec2 v2(1.1, 1.2);
  Vec2 v3 = v1 - v2;
  EXPECT_NEAR(v3[0], 0.0, EPS);
  EXPECT_NEAR(v3[1], 1.0, EPS);
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
TEST(Vec3Test, Minus) {
  Vec3 v1(1.1, 2.2, 3.3);
  Vec3 v2(1.1, 1.2, 3.0);
  Vec3 v3 = v1 - v2;
  EXPECT_NEAR(v3[0], 0.0, EPS);
  EXPECT_NEAR(v3[1], 1.0, EPS);
  EXPECT_NEAR(v3[2], 0.3, EPS);
}
TEST(Vec3Test, Multipy) {
  Vec3 v1(1.0, 2.2, 3.3);
  Vec3 v2 = v1 * 1.1;
  EXPECT_NEAR(v2[0], 1.1, EPS);
  EXPECT_NEAR(v2[1], 2.2 * 1.1, EPS);
  EXPECT_NEAR(v2[2], 3.3 * 1.1, EPS);
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
TEST(VoxelTest, IntersectsSimple) {
  Voxel v({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  Ray r({2.0, 2.0, 2.0}, {1.0, 1.0, 1.0});
  EXPECT_FALSE(v.intersects(r));
  r = Ray({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0});
  EXPECT_TRUE(v.intersects(r));
  r = Ray({2.0, 2.0, 2.0}, {-1.0, -1.0, -1.0});
  EXPECT_TRUE(v.intersects(r));
}
Vec3 randomVec3(double lo, double hi, std::mt19937& random_engine) {
  assert(lo < hi);
  std::uniform_real_distribution dist(lo, hi);
  return Vec3(dist(random_engine), dist(random_engine), dist(random_engine));
}
double randomLogUniformReal(double lo_log, double hi_log,
                            std::mt19937& random_engine) {
  assert(lo_log < hi_log);
  std::uniform_real_distribution dist(lo_log, hi_log);
  return std::pow(2.0, dist(random_engine));
}
// Random test for rays intersecting inner points of the voxel
TEST(VoxelTest, IntersectsRandom) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    std::uniform_real_distribution<double> dist(-1, 1);
    Voxel v(randomVec3(-1, 1, mt), randomVec3(-1, 1, mt));

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
      Vec3 other_point = randomVec3(-1, 1, mt);
      double scale = randomLogUniformReal(-10, 10, mt);
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
  v.cover(&t);
  Voxel c({-2.2, -1.0, -1.2}, {5.2, 0.0, 1.0});
  EXPECT_NEAR((v.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((v.hi - c.hi).norm(), 0, EPS);
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
Vec2 randomBaryCoords(std::mt19937& random_engine) {
  std::uniform_real_distribution<double> dist(0, 1);
  Vec2 v(dist(random_engine), dist(random_engine));
  // limit sum to 1
  v[1] = std::min(1 - v[0], v[1]);
  return v;
}
// tests that random barycentric coordinates will be on the triangle plane
TEST(TriangleTest, PointFromBaryRandom) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    Triangle t(randomVec3(-1.0, 1.0, mt), randomVec3(-1.0, 1.0, mt),
               randomVec3(-1.0, 1.0, mt));
    Vec2 bary_coords(randomBaryCoords(mt));
    Vec3 p(t.pointFromBary(bary_coords));
    Vec3 dot_normal = (p - t.p0).cross(t.p1 - t.p0).dot(t.p2 - t.p0);
    ASSERT_NEAR(dot_normal.norm(), 0, EPS);
  }
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
TEST(TriangleTest, RayIntersectionRandom) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 100'000; ++i) {
    Triangle random_t(randomVec3(-1.0, 1.0, mt), randomVec3(-1.0, 1.0, mt),
                      randomVec3(-1.0, 1.0, mt));

    Vec2 bary_coords = randomBaryCoords(mt);
    // skip coordinates that are too close to border
    if (bary_coords[0] < EPS || bary_coords[1] < EPS ||
        bary_coords[0] + bary_coords[1] > 1 - EPS) {
      continue;
    }
    Vec3 inner_point(random_t.pointFromBary(bary_coords));
    Vec3 ray_origin(randomVec3(-1.0, 1.0, mt));
    // ray origin should not be on the plane
    // TODO implement this to geometry.cpp?
    Vec3 tmp = (ray_origin - random_t.p0).cross(random_t.p1 - random_t.p0);
    if (std::abs(tmp.dot(random_t.p2 - random_t.p0)) < EPS) {
      continue;
    }
    double scale = randomLogUniformReal(-1.0, 1.0, mt);
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
TEST(PlanePolygon, ConstructorPoints) {
  PlanePolygon p(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  std::vector<Vec3> r = {p[0], p[1], p[2]};
  EXPECT_NEAR((p[0] - Vec3(1.0, 0.0, 0.0)).norm(), 0, EPS);
  EXPECT_NEAR((p[1] - Vec3(1.1)).norm(), 0, EPS);
  EXPECT_NEAR((p[2] - Vec3(2.2)).norm(), 0, EPS);
}
TEST(PlanePolygon, ConstructorTriangle) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(&t);
  std::vector<Vec3> r = {p[0], p[1], p[2]};
  EXPECT_NEAR((p[0] - Vec3(1.0, 0.0, 0.0)).norm(), 0, EPS);
  EXPECT_NEAR((p[1] - Vec3(1.1)).norm(), 0, EPS);
  EXPECT_NEAR((p[2] - Vec3(2.2)).norm(), 0, EPS);
}
TEST(PlanePolygon, BoundingBoxSimple) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(&t);
  Voxel r(p.getBoundingBox());
  Voxel c(Vec3(1.0, 0.0, 0.0), Vec3(2.2));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  p = PlanePolygon(&t);
  r = Voxel(p.getBoundingBox());
  c = Voxel(Vec3(-12.23, -15.2, 2.2), Vec3(1.1, 2.2, 3.9999));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);
}
TEST(PlanePolygon, IntersectSimple) {
  Triangle t(Vec3(0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 1.0));
  PlanePolygon p(&t);
  p.intersect({2, 0.5}, 0);
  Voxel r = p.getBoundingBox();
  Voxel c = Voxel(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.5, 0.5));
  EXPECT_NEAR((r.lo - c.lo).norm(), 0, EPS);
  EXPECT_NEAR((r.hi - c.hi).norm(), 0, EPS);
}
// TODO add more intersection tests
TEST(PlanePolygon, Size) {
  Triangle t(Vec3(1.0, 0.0, 0.0), Vec3(1.1), Vec3(2.2));
  PlanePolygon p(&t);
  EXPECT_EQ(p.size(), 3);
}
TEST(ClipTriangle, MinMin) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(&t);
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
  ClipTriangle ct(&t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  t = Triangle(Vec3(1.0, 0.0, 2.5), Vec3(1.0, 2.2, 3.9999),
               Vec3(1.0, -15.2, 2.2));
  ct = ClipTriangle(&t);
  EXPECT_TRUE(ct.isAxisAligned(0));
  EXPECT_FALSE(ct.isAxisAligned(1));
  EXPECT_FALSE(ct.isAxisAligned(2));

  t = Triangle(Vec3(1.0, 2.0, 12), Vec3(2.0, 2.0, 3.9999), Vec3(3.0, 2.0, 2.2));
  ct = ClipTriangle(&t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  EXPECT_TRUE(ct.isAxisAligned(1));
  EXPECT_FALSE(ct.isAxisAligned(2));

  t = Triangle(Vec3(1.0, 12, 3.3), Vec3(2.0, 2.0, 3.3), Vec3(3.0, 3.0, 3.3));
  ct = ClipTriangle(&t);
  EXPECT_FALSE(ct.isAxisAligned(0));
  EXPECT_FALSE(ct.isAxisAligned(1));
  EXPECT_TRUE(ct.isAxisAligned(2));
}
TEST(ClipTriangle, OverlapSidesSimple) {
  Triangle t(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
             Vec3(-12.23, -15.2, 2.2));
  ClipTriangle ct(&t);
  std::pair<bool, bool> r = ct.overlapsSides({1, 1.0}, 1);
  EXPECT_EQ(r, std::make_pair(true, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(&t);
  r = ct.overlapsSides({0, -100.0}, 1);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(-0.1, 0.0, 2.5), Vec3(1.1, 2.2, 3.9999),
               Vec3(-12.23, -15.2, 2.2));
  ct = ClipTriangle(&t);
  r = ct.overlapsSides({2, 100.0}, 0);
  EXPECT_EQ(r, std::make_pair(true, false));
}
TEST(ClipTriangle, OverlapSidesParallel) {
  Triangle t(Vec3(0.0, 0.0, 2.5), Vec3(0.0, 2.2, 3.9999),
             Vec3(0.0, -15.2, 2.2));
  ClipTriangle ct(&t);
  std::pair<bool, bool> r = ct.overlapsSides({0, 0.0}, 1);
  EXPECT_EQ(r, std::make_pair(false, true));

  t = Triangle(Vec3(1.0, 1.1, 2.5), Vec3(1.2, 1.1, 3.9999),
               Vec3(1.1, 1.1, 2.2));
  ct = ClipTriangle(&t);
  r = ct.overlapsSides({1, 1.1}, 0);
  EXPECT_EQ(r, std::make_pair(true, false));
}
TEST(ClipTriangle, ClipSimple) {
  Triangle t(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ClipTriangle ct(&t);
  ct.clip({0, 0.5}, 0);
  EXPECT_NEAR(ct.max(0), 0.5, EPS);
  EXPECT_NEAR(ct.max(1), 1.0, EPS);
  EXPECT_NEAR(ct.min(0), 0.0, EPS);

  t = Triangle(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ct = ClipTriangle(&t);
  ct.clip({0, 0.5}, 1);
  EXPECT_NEAR(ct.min(0), 0.5, EPS);
  EXPECT_NEAR(ct.max(1), 0.5, EPS);
}
TEST(ClipTriangle, ClipParallel) {
  Triangle t(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
  ClipTriangle ct(&t);
  ct.clip({2, 0.0}, 0);
  EXPECT_NEAR(ct.min(2), 0.0, EPS);
  EXPECT_NEAR(ct.max(2), 0.0, EPS);
  EXPECT_NEAR(ct.min(0), 0.0, EPS);
  EXPECT_NEAR(ct.max(0), 1.0, EPS);
}
}  // namespace
