#include "../../kd_tree.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <chrono>
#include <tuple>

#include "../../geometry.h"
#include "../../raycaster.h"
#include "../test_utils.h"

using ::testing::ContainerEq;
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
MATCHER_P(TriangleVecEq, v, "should equal " + PrintToString(v)) {
  if (v.size() != arg.size()) return 0;
  for (size_t i = 0; i < v.size(); ++i) {
    if ((v[i].p0 - arg[i].p0).norm() > EPS) return 0;
    if ((v[i].p1 - arg[i].p1).norm() > EPS) return 0;
    if ((v[i].p2 - arg[i].p2).norm() > EPS) return 0;
  }
  return 1;
}
TEST(RelativeSubvoxelAreas, Simple) {
  Voxel v(Vec3(0.0), Vec3(1.0));
  AxisPlane ap = {0, 0.5};
  double l_area, r_area;
  std::tie(l_area, r_area) = relativeSubvoxelAreas(v, ap);
  double c = (2 + 0.5 * 4) / 6.0;
  EXPECT_NEAR(l_area, c, EPS);
  EXPECT_NEAR(r_area, c, EPS);

  v = Voxel(Vec3(0.0), Vec3(1.0));
  ap = {0, 0.0};
  std::tie(l_area, r_area) = relativeSubvoxelAreas(v, ap);
  double l_c = 2.0 / 6.0;
  double r_c = 1.0;
  EXPECT_NEAR(l_area, l_c, EPS);
  EXPECT_NEAR(r_area, r_c, EPS);
}
TEST(RelativeSubvoxelAreas, FlatVoxel) {
  Voxel v(Vec3(0.0), Vec3(1.0, 1.0, 0.0));
  AxisPlane ap = {2, 0};
  double l_area, r_area;
  std::tie(l_area, r_area) = relativeSubvoxelAreas(v, ap);
  EXPECT_NEAR(l_area, 1.0, EPS);
  EXPECT_NEAR(r_area, 1.0, EPS);
}
TEST(SurfaceAreaHeuristic, Simple) {
  double l_area = 0.5;
  double r_area = 0.5;
  TriangleCounts counts;
  counts.left = 10;
  counts.plane = 0;
  counts.right = 2;
  double traversal_cost = 10.0;
  double intersection_cost = 1.0;
  double cost;
  bool side;
  std::tie(cost, side) =
      surfaceAreaHeuristic(l_area, r_area, counts,
                           traversal_cost, intersection_cost);
  double correct =
      traversal_cost + intersection_cost * (l_area * counts.left + r_area * counts.right);
  EXPECT_NEAR(cost, correct, EPS);
}
TEST(SurfaceAreaHeuristic, TrianglesOnPlane) {
  double l_area = 0.5;
  double r_area = 0.5;
  TriangleCounts counts;
  counts.left = 10;
  counts.plane = 4;
  counts.right = 2;
  double traversal_cost = 10.0;
  double intersection_cost = 1.0;
  double cost;
  bool side;

  std::tie(cost, side) =
      surfaceAreaHeuristic(l_area, r_area, counts,
                           traversal_cost, intersection_cost);
  double l_correct =
      traversal_cost +
      intersection_cost * (l_area * (counts.left + counts.plane) + r_area * counts.right);
  double r_correct =
      traversal_cost +
      intersection_cost * (l_area * counts.left + r_area * (counts.right + counts.plane));

  EXPECT_NEAR(cost, std::min(l_correct, r_correct), EPS);
  EXPECT_EQ(side, l_correct <= r_correct ? 0 : 1);
}
/* Disable at least for time being

// tests that the extractTriangles reverses createClipTriangles correctly
TEST(CreateAndExtractTriangles, Simple) {
  std::mt19937 mt(1337);
  Vector<Triangle> tv = test::randomTriangleVector(-10, 10, 10, 100, mt);
  Vector<ClipTriangle> ct = createClipTriangles(tv);
  Vector<Triangle> restored = extractTriangles(ct);
  EXPECT_THAT(restored, TriangleVecEq(tv));
}
*/
// Test that tree building and queries work for an empty tree
TEST(TreeBuilder, Empty) {
  Vector<SceneTriangle*> scene_p;
  Tree t = buildKdTree(scene_p, 1.0, 40.0);
  Ray r(Vec3(1.0), Vec3(1.0));
  ScenePoint sp = t.getClosestRayIntersection(r);
  ASSERT_EQ(sp.scene_triangle, nullptr);
}
// tests both the TreeBuilder and Tree by building a tree and
// making queries
TEST(TreeBuilderKdTreeQueries, Random3d) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 1000; ++i) {
    double triangle_scale = test::randomLogUniformReal(-4, 7, mt);
    Vector<SceneTriangle> scene =
        test::randomSceneTriangleVector(0.0, 100.0, triangle_scale, 100, mt);
    Vector<SceneTriangle*> scene_p;
    for (auto& x : scene) {
      scene_p.pushBack(&x);
    }
    Vector<Triangle> triangles = extractTriangles(scene_p);
    Tree t = buildKdTree(scene_p, 1.0, 40.0);
    for (int j = 0; j < 100; ++j) {
      Ray r(test::randomVec3(0.0, 100.0, mt), test::randomVec3(-1.0, 1.0, mt));
      ScenePoint sp = t.getClosestRayIntersection(r);
      RayTriangleIntersection correct =
          firstRayTriangleIntersection(triangles, r);
      SceneTriangle* correct_st = nullptr;
      if (correct.index < scene_p.size()) correct_st = scene_p[correct.index];
      ASSERT_EQ(sp.scene_triangle, correct_st);
      ASSERT_THAT(sp.bary_coords, VecEq(correct.bary_coords));
    }
  }
}
TEST(TreeBuilderKdTreeQueries, Random3dAxisParallel) {
  std::mt19937 mt(1337);
  for (int i = 0; i < 1000; ++i) {
    double triangle_scale = test::randomLogUniformReal(-4, 7, mt);
    Vector<SceneTriangle> scene =
        test::randomSceneTriangleVector(0.0, 100.0, triangle_scale, 100, mt);
    // make triangles parallel to some axis
    std::uniform_int_distribution axis_dist(0, 2);
    for (size_t j = 0; j < scene.size(); ++j) {
      int ax = axis_dist(mt);
      scene[j].triangle.p1[ax] = scene[j].triangle.p0[ax];
      scene[j].triangle.p2[ax] = scene[j].triangle.p0[ax];
    }
    Vector<SceneTriangle*> scene_p;
    for (auto& x : scene) {
      scene_p.pushBack(&x);
    }
    Vector<Triangle> triangles = extractTriangles(scene_p);
    Tree t = buildKdTree(scene_p, 1.0, 40.0);
    for (int j = 0; j < 100; ++j) {
      Ray r(test::randomVec3(0.0, 100.0, mt), test::randomVec3(-1.0, 1.0, mt));
      ScenePoint sp = t.getClosestRayIntersection(r);
      RayTriangleIntersection correct =
          firstRayTriangleIntersection(triangles, r);
      SceneTriangle* correct_st = nullptr;
      if (correct.index < scene_p.size()) correct_st = scene_p[correct.index];
      ASSERT_EQ(sp.scene_triangle, correct_st);
      ASSERT_THAT(sp.bary_coords, VecEq(correct.bary_coords));
    }
  }
}
// tests the TreeBuilder by checking that the resulting kd tree
// has good enough performance
// (used as a proxy for a correctly built kd tree)
TEST(TreeBuilderKdTreeQueries, Random3dSmallTrianglePerformance) {
  std::mt19937 mt(1337);
  std::chrono::duration<double, std::milli> kd_tree_time(0.0);
  std::chrono::duration<double, std::milli> naive_time(0.0);
  int n_rays = 0;
  int n_intersections = 0;
  for (int i = 0; i < 5; ++i) {
    Vector<SceneTriangle> scene =
        test::randomSceneTriangleVector(0.0, 100.0, 0.05, 100000, mt);
    Vector<SceneTriangle*> scene_p;
    for (auto& x : scene) {
      scene_p.pushBack(&x);
    }
    Vector<Triangle> triangles = extractTriangles(scene_p);
    Tree t = buildKdTree(scene_p, 1.0, 40.0);
    std::cerr << "tree " << i << " done " << std::endl;
    for (int i = 0; i < 100; ++i) {
      Ray r(test::randomVec3(0.0, 1.0, mt), test::randomVec3(0.0, 1.0, mt));
      auto t0 = std::chrono::high_resolution_clock::now();
      ScenePoint sp = t.getClosestRayIntersection(r);
      auto t1 = std::chrono::high_resolution_clock::now();
      RayTriangleIntersection correct =
          firstRayTriangleIntersection(triangles, r);
      auto t2 = std::chrono::high_resolution_clock::now();
      kd_tree_time += t1 - t0;
      naive_time += t2 - t1;
      ++n_rays;
      SceneTriangle* correct_st = nullptr;
      if (correct.index < scene_p.size()) correct_st = scene_p[correct.index];
      ASSERT_EQ(sp.scene_triangle, correct_st);
      ASSERT_THAT(sp.bary_coords, VecEq(correct.bary_coords));
    }
  }
  EXPECT_LE(30 * kd_tree_time, naive_time);
  std::cerr << "kd_tree_time: " << kd_tree_time.count() << std::endl;
  std::cerr << "naive_time: " << naive_time.count() << std::endl;
  std::cerr << "total rays: " << n_rays << std::endl;
  std::cerr << "intersections: " << n_intersections << std::endl;
}
TEST(TreeBuilderKdTreeQueries, Random3dLargeTrianglePerformance) {
  std::mt19937 mt(1337);
  std::chrono::duration<double, std::milli> kd_tree_time(0.0);
  std::chrono::duration<double, std::milli> naive_time(0.0);
  int n_rays = 0;
  int n_intersections = 0;
  for (int i = 0; i < 5; ++i) {
    Vector<SceneTriangle> scene =
        test::randomSceneTriangleVector(0.0, 100.0, 3, 20000, mt);
    Vector<SceneTriangle*> scene_p;
    for (auto& x : scene) {
      scene_p.pushBack(&x);
    }
    Vector<Triangle> triangles = extractTriangles(scene_p);
    Tree t = buildKdTree(scene_p, 1.0, 40.0);
    std::cerr << "tree " << i << " done " << std::endl;
    for (int i = 0; i < 100; ++i) {
      Ray r(test::randomVec3(0.0, 1.0, mt), test::randomVec3(0.0, 1.0, mt));
      auto t0 = std::chrono::high_resolution_clock::now();
      ScenePoint sp = t.getClosestRayIntersection(r);
      auto t1 = std::chrono::high_resolution_clock::now();
      RayTriangleIntersection correct =
          firstRayTriangleIntersection(triangles, r);
      auto t2 = std::chrono::high_resolution_clock::now();
      kd_tree_time += t1 - t0;
      naive_time += t2 - t1;
      ++n_rays;
      SceneTriangle* correct_st = nullptr;
      if (correct.index < scene_p.size()) correct_st = scene_p[correct.index];
      ASSERT_EQ(sp.scene_triangle, correct_st);
      ASSERT_THAT(sp.bary_coords, VecEq(correct.bary_coords));
    }
  }
  EXPECT_LE(30 * kd_tree_time, naive_time);
  std::cerr << "kd_tree_time: " << kd_tree_time.count() << std::endl;
  std::cerr << "naive_time: " << naive_time.count() << std::endl;
  std::cerr << "total rays: " << n_rays << std::endl;
  std::cerr << "intersections: " << n_intersections << std::endl;
}
TEST(TreeBuilderKdTreeQueries, TestTrianglesIntersectSegmentSimple) {
  Triangle triangle(Vec3(0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0));
  Vector<SceneTriangle> scene;
  scene.pushBack({triangle, {}, {}});
  Vector<SceneTriangle*> scene_p;
  scene_p.pushBack(&scene[0]);
  Tree tree = buildKdTree(scene_p, 1.0, 5.0);
  Vec3 a(-0.5, 0.2, 0.2);
  Vec3 b(0.5, 0.2, 0.2);
  EXPECT_TRUE(tree.trianglesIntersectSegment(a, b));
  a = Vec3(0.1, 0.2, 0.2);
  b = Vec3(0.5, 0.2, 0.2);
  EXPECT_FALSE(tree.trianglesIntersectSegment(a, b));
  a = Vec3(0.1, 1.2, 1.2);
  b = Vec3(0.5, 1.2, 1.2);
  EXPECT_FALSE(tree.trianglesIntersectSegment(a, b));
}

}  // namespace
