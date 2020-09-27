#include <random>

#include "geometry.h"
#include "raycaster.h"

const double TEST_EPS = 1e-7;

namespace test {

Vec3 randomVec3(double lo, double hi, std::mt19937& random_engine);
double randomLogUniformReal(double lo_log, double hi_log,
                            std::mt19937& random_engine);
Vec2 randomBaryCoords(std::mt19937& random_engine);
std::vector<Triangle> randomTriangleVector(double lo, double hi,
                                           double max_triangle_size, int n,
                                           std::mt19937 random_engine);
std::vector<SceneTriangle> randomSceneTriangleVector(
    double lo, double hi, double max_triangle_size, int n,
    std::mt19937 random_engine);
bool pointOnTrianglePlane(const Triangle& t, const Vec3& p);

std::vector<int> randomIntVector(int lo, int hi, int n,
                                 std::mt19937& random_engine);

}  // namespace test
