#include "test_utils.h"
namespace test {
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
Vec2 randomBaryCoords(std::mt19937& random_engine) {
  std::uniform_real_distribution<double> dist(0, 1);
  Vec2 v(dist(random_engine), dist(random_engine));
  // limit sum to 1
  v[1] = std::min(1 - v[0], v[1]);
  return v;
}
// generates triangles into area [lo-max_triangle_size, hi+max_triangle_size]
std::vector<Triangle*> randomTriangleVector(double lo, double hi,
                                            double max_triangle_size, int n,
                                            std::mt19937 random_engine) {
  assert(lo < hi);
  std::vector<Triangle*> v;
  for (int i = 0; i < n; ++i) {
    Vec3 p0(randomVec3(lo, hi, random_engine));
    Vec3 p1(p0 +
            randomVec3(-max_triangle_size, max_triangle_size, random_engine));
    Vec3 p2(p0 +
            randomVec3(-max_triangle_size, max_triangle_size, random_engine));
    Triangle* t = new Triangle(p0, p1, p2);
    // check if the triangle has some area
    if ((t->p1 - t->p0).cross(t->p2 - t->p0).norm() > EPS) {
      v.push_back(t);
    }
  }
  return v;
}
bool pointOnTrianglePlane(const Triangle& t, const Vec3& p) {
  return std::abs((p - t.p0).cross(t.p1 - t.p0).dot(t.p2 - t.p0)) < TEST_EPS;
}

}  // namespace test
