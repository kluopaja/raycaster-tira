#include "kd_tree.h"
#include "test_utils.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <chrono>
#include <tuple>

using ::testing::ContainerEq;
using ::testing::PrintToString;
namespace {
    TEST(RelativeSubvoxelAreas, Simple) {
        Voxel v(Vec3(0.0), Vec3(1.0));
        AxisPlane ap = {0, 0.5};
        double l_area, r_area;
        std::tie(l_area, r_area) = relativeSubvoxelAreas(v, ap);
        double c = (2 + 0.5*4)/6.0;
        EXPECT_NEAR(l_area, c, EPS);
        EXPECT_NEAR(r_area, c, EPS);

        v = Voxel(Vec3(0.0), Vec3(1.0));
        ap = {0, 0.0};
        std::tie(l_area, r_area) = relativeSubvoxelAreas(v, ap);
        double l_c = 2.0/6.0;
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
        int n_left = 10;
        int n_plane = 0;
        int n_right = 2;
        double traversal_cost = 10.0;
        double intersection_cost = 1.0;
        double cost;
        bool side;
        std::tie(cost, side)  = surfaceAreaHeuristic(l_area, r_area, n_left,
                                                  n_plane, n_right,
                                                  traversal_cost,
                                                  intersection_cost);
        double correct = traversal_cost
                         + intersection_cost
                           * (l_area * n_left + r_area * n_right);
        EXPECT_NEAR(cost, correct, EPS);
    }
    TEST(SurfaceAreaHeuristic, TrianglesOnPlane) {
        double l_area = 0.5;
        double r_area = 0.5;
        int n_left = 10;
        int n_plane = 4;
        int n_right = 2;
        double traversal_cost = 10.0;
        double intersection_cost = 1.0;
        double cost;
        bool side;
        std::tie(cost, side)  = surfaceAreaHeuristic(l_area, r_area, n_left,
                                                  n_plane, n_right,
                                                  traversal_cost,
                                                  intersection_cost);
        double l_correct = traversal_cost
                         + intersection_cost
                           * (l_area * (n_left + n_plane) + r_area * n_right);
        double r_correct = traversal_cost
                         + intersection_cost
                           * (l_area * n_left + r_area * (n_right + n_plane));

        EXPECT_NEAR(cost, std::min(l_correct, r_correct), EPS);
        EXPECT_EQ(side, l_correct <= r_correct ? 0 : 1);
// tests both the TreeBuilder and Tree by building a tree and
// making queries
TEST(TreeBuilderKdTreeQueries, Random3d) {
  std::mt19937 mt(1337);
  int n_rays = 0;
  int n_same_answer = 0;
  for (int i = 0; i < 1000; ++i) {
    double triangle_scale = test::randomLogUniformReal(-4, 7, mt);
    std::vector<Triangle*> scene =
        test::randomTriangleVector(0.0, 100.0, triangle_scale, 100, mt);
    Tree t = buildKdTree(scene, 1.0, 40.0);
    for (int j = 0; j < 100; ++j) {
      Ray r(test::randomVec3(0.0, 1.0, mt), test::randomVec3(0.0, 1.0, mt));
      TrianglePoint tp = t.getClosestRayIntersection(r);
      TrianglePoint correct = firstRayTriangleIntersection(scene, r);
      ++n_rays;
      if (tp.triangle == correct.triangle) ++n_same_answer;
      EXPECT_THAT(tp.bary_coords, VecEq(correct.bary_coords));
// tests the TreeBuilder by checking that the resulting kd tree
// has good enough performance
// (used as a proxy for a correctly built kd tree)
TEST(TreeBuilderKdTreeQueries, Random3dPerformance) {
  std::mt19937 mt(1337);
  std::chrono::duration<double, std::milli> kd_tree_time(0.0);
  std::chrono::duration<double, std::milli> naive_time(0.0);
  int n_rays = 0;
  int n_same_answer = 0;
  for (int i = 0; i < 5; ++i) {
    std::vector<Triangle*> scene =
        test::randomTriangleVector(0.0, 100.0, 0.05, 100000, mt);
    Tree t = buildKdTree(scene, 1.0, 40.0);
    std::cerr << "tree " << i << " done " << std::endl;
    for (int i = 0; i < 100; ++i) {
      Ray r(test::randomVec3(0.0, 1.0, mt), test::randomVec3(0.0, 1.0, mt));
      auto t0 = std::chrono::high_resolution_clock::now();
      TrianglePoint tp = t.getClosestRayIntersection(r);
      auto t1 = std::chrono::high_resolution_clock::now();
      TrianglePoint correct = firstRayTriangleIntersection(scene, r);
      auto t2 = std::chrono::high_resolution_clock::now();
      kd_tree_time += t1 - t0;
      naive_time += t2 - t1;
      ++n_rays;
      if (tp.triangle == correct.triangle) ++n_same_answer;
      ASSERT_EQ(tp.triangle, correct.triangle);
      EXPECT_THAT(tp.bary_coords, VecEq(correct.bary_coords));
    }
    // tests that the extractTriangles reverses createClipTriangles correctly
    TEST(CreateAndExtractTriangles, Simple) {
        std::mt19937 mt(1337);
        std::vector<Triangle*> tv = randomTriangleVector(-10, 10, 100, mt);
        std::vector<ClipTriangle> ct = createClipTriangles(tv);
        std::vector<Triangle*> restored = extractTriangles(ct);
        EXPECT_THAT(restored, ContainerEq(tv));
    }
    /*
 class Scene2D: public ::testing::Test {
     protected:
         std::vector<Triangle*> scene;
         void SetUp() override {
             scene.push_back(new Triangle(Vec3(0.0, 0.0, 0.0),
                                          Vec3(1.0, 0.0, 0.0),
                                          Vec3(0.0, 1.0, 0.0)));
             scene.push_back(new Triangle(Vec3(2.0, 2.0, 0.0),
                                          Vec3(2.5, 2.0, 0.0),
                                          Vec3(2.0, 2.5, 0.0)));
             scene.push_back(new Triangle(Vec3(2.2, 1.0, 0.0),
                                          Vec3(1.3, 2.2, 0.0),
                                          Vec3(1.1, 2.1, 0.0)));
         }
         void TearDown() override {
             for(int i = 0; i < scene.size(); ++i) {
                 delete scene[i];
             }
         }
};
    TEST_F(Scene2D, BuildKdTree) {
        Tree t = buildKdTree(scene, 1.0, 2.0);
    }
    */
} // namespace
  }
  EXPECT_LE(30 * kd_tree_time, naive_time);
  std::cerr << "kd_tree_time: " << kd_tree_time.count() << std::endl;
  std::cerr << "naive_time: " << naive_time.count() << std::endl;
  std::cerr << "total rays: " << n_rays << std::endl;
}

}  // namespace
