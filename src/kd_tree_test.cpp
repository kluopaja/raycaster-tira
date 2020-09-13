#include "kd_tree.h"
#include "test_utils.h"

#include <tuple>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using ::testing::ContainerEq;
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
