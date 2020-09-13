#include "test_utils.h"
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
std::vector<Triangle*> randomTriangleVector(double lo, double hi, int n,
                                            std::mt19937 random_engine) {
    assert(lo < hi);
    std::vector<Triangle*> v;
    for(int i = 0; i < n; ++i) {
        Triangle* t = new Triangle(randomVec3(lo, hi, random_engine),
                                   randomVec3(lo, hi, random_engine),
                                   randomVec3(lo, hi, random_engine));
        // check if the triangle has some area
        if((t->p1 - t->p0).cross(t->p2 - t->p0).norm() > EPS) {
            v.push_back(t);
        }
    }
    return v;
}
