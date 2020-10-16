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
Vector<Triangle> randomTriangleVector(double lo, double hi,
                                           double max_triangle_size, int n,
                                           std::mt19937 random_engine) {
  assert(lo < hi);
  Vector<Triangle> v;
  for (int i = 0; i < n; ++i) {
    Vec3 p0(randomVec3(lo, hi, random_engine));
    Vec3 p1(p0 +
            randomVec3(-max_triangle_size, max_triangle_size, random_engine));
    Vec3 p2(p0 +
            randomVec3(-max_triangle_size, max_triangle_size, random_engine));
    Triangle t(p0, p1, p2);
    // check if the triangle has some area
    if ((t.p1 - t.p0).cross(t.p2 - t.p0).norm() > EPS) {
      v.pushBack(t);
    }
  }
  return v;
}
// generates scene triangles into area [lo-max_triangle_size,
// hi+max_triangle_size]
Vector<SceneTriangle> randomSceneTriangleVector(
    double lo, double hi, double max_triangle_size, int n,
    std::mt19937 random_engine) {
  assert(lo < hi);
  Vector<Triangle> t =
      randomTriangleVector(lo, hi, max_triangle_size, n, random_engine);
  Vector<SceneTriangle> v;
  for (size_t i = 0; i < t.size(); ++i) {
    SceneTriangle scene_triangle;
    scene_triangle.triangle = t[i];
    v.pushBack(scene_triangle);
  }
  return v;
}
bool pointOnTrianglePlane(const Triangle& t, const Vec3& p) {
  return std::abs((p - t.p0).cross(t.p1 - t.p0).dot(t.p2 - t.p0)) < TEST_EPS;
}

Vector<int> randomIntVector(int lo, int hi, int n,
                                 std::mt19937& random_engine) {
  assert(lo <= hi);
  assert(n >= 0);
  std::uniform_int_distribution dist(lo, hi);
  Vector<int> out;
  for (int i = 0; i < n; ++i) {
    out.pushBack(dist(random_engine));
  }
  return out;
}

Material randomMaterial(std::mt19937& random_engine) {
  Material m = Material();
  m.emitted = randomVec3(0.01, 1.0, random_engine);
  m.diffuse = randomVec3(0.01, 1.0, random_engine);
  m.specular = randomVec3(0.01, 1.0, random_engine);
  // cap to 1
  for (int i = 0; i < 3; ++i) {
    if(m.diffuse[i] + m.specular[i] > 1) {
      m.specular[i] = 1 - m.diffuse[i];
    }
  }
  m.transparent = randomVec3(0.01, 1.0, random_engine);
  std::uniform_real_distribution dist(0.1, 100.0);
  m.specular_exp = dist(random_engine);
  m.index_of_refraction = dist(random_engine);
  return m;
}
double mean(const Vector<double>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double variance(const Vector<double>& v) {
  double mean = test::mean(v);
  double sum_of_squares = 0;
  for(double x: v) {
    sum_of_squares += std::pow(x - mean, 2);
  }
  return sum_of_squares / v.size();
}
double standard_error_mean(const Vector<double>& v) {
  return std::sqrt(variance(v) / v.size());
}
}  // namespace test
