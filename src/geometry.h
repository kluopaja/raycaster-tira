#ifndef RAYCASTER_GEOMETRY_H
#define RAYCASTER_GEOMETRY_H
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

const double INF = std::numeric_limits<double>::infinity();
const double EPS = 1e-10;
const double kPi = std::atan(1) * 4;

class Triangle;
class Ray;
class AxisPlane;

// Implement Vec2 and Vec3 separately at least for now
class Vec2 {
 public:
  Vec2() = default;
  Vec2(double x);
  Vec2(double x, double y);
  double dot(const Vec2& b) const;
  double sum() const;
  // element-wise product
  Vec2 multiply(const Vec2& b) const;
  friend Vec2& operator+=(Vec2& a, const Vec2& b);
  friend Vec2 operator+(const Vec2& a, const Vec2& b);
  friend Vec2& operator-=(Vec2& a, const Vec2& b);
  friend Vec2 operator-(const Vec2& a, const Vec2& b);
  friend Vec2& operator*=(Vec2& a, const double b);
  friend Vec2 operator*(const Vec2& a, const double b);
  friend Vec2 operator*(const double b, const Vec2& a);
  friend Vec2& operator/=(Vec2& a, const double b);
  friend Vec2 operator/(const Vec2& a, const double b);
  double& operator[](int index);
  double operator[](int index) const;
  double norm() const;
  friend std::ostream& operator<<(std::ostream& out, const Vec2& a);

 private:
  double v[2];
};
inline double Vec2::dot(const Vec2& b) const {
  return v[0] * b.v[0] + v[1] * b.v[1];
}
inline Vec2 Vec2::multiply(const Vec2& b) const {
  return Vec2(v[0] * b.v[0], v[1] * b.v[1]);
}
inline double Vec2::sum() const {
  return v[0] + v[1];
}
inline Vec2& operator+=(Vec2& a, const Vec2& b) {
  a.v[0] += b.v[0];
  a.v[1] += b.v[1];
  return a;
}
inline Vec2 operator+(const Vec2& a, const Vec2& b) {
  return Vec2(a.v[0] + b.v[0], a.v[1] + b.v[1]);
}
inline Vec2& operator-=(Vec2& a, const Vec2& b) {
  a.v[0] -= b.v[0];
  a.v[1] -= b.v[1];
  return a;
}
inline Vec2 operator-(const Vec2& a, const Vec2& b) {
  return Vec2(a.v[0] - b.v[0], a.v[1] - b.v[1]);
}
inline Vec2& operator*=(Vec2& a, double b) {
  a.v[0] *= b;
  a.v[1] *= b;
  return a;
}
inline Vec2 operator*(const Vec2& a, const double b) {
  return Vec2(b * a.v[0], b * a.v[1]);
}
inline Vec2 operator*(const double b, const Vec2& a) {
  return Vec2(b * a.v[0], b * a.v[1]);
}
inline Vec2& operator/=(Vec2& a, double b) {
  a.v[0] /= b;
  a.v[1] /= b;
  return a;
}
inline Vec2 operator/(const Vec2& a, const double b) {
  return Vec2(a.v[0] / b, a.v[1] / b);
}
inline double& Vec2::operator[](int index) { return v[index]; }
inline double Vec2::operator[](int index) const { return v[index]; }
inline double Vec2::norm() const { return std::sqrt(this->dot(*this)); }
class Vec3 {
 public:
  Vec3() = default;
  Vec3(double x);
  Vec3(double x, double y, double z);
  double dot(const Vec3& b) const;
  Vec3 cross(const Vec3& b) const;
  double sum() const;
  // element-wise product
  Vec3 multiply(const Vec3& b) const;
  friend Vec3& operator+=(Vec3& a, const Vec3& b);
  friend Vec3 operator+(const Vec3& a, const Vec3& b);
  friend Vec3& operator-=(Vec3& a, const Vec3& b);
  friend Vec3 operator-(const Vec3& a, const Vec3& b);
  friend Vec3& operator*=(Vec3& a, const double b);
  friend Vec3 operator*(const Vec3& a, const double b);
  friend Vec3 operator*(const double b, const Vec3& a);
  friend Vec3& operator/=(Vec3& a, const double b);
  friend Vec3 operator/(const Vec3& a, const double b);
  double& operator[](int index);
  double operator[](int index) const;
  double norm() const;
  friend std::ostream& operator<<(std::ostream& out, const Vec3& a);

 private:
  double v[3];
};
inline double Vec3::dot(const Vec3& b) const {
  return v[0] * b.v[0] + v[1] * b.v[1] + v[2] * b.v[2];
}
inline Vec3 Vec3::cross(const Vec3& b) const {
  return {v[1] * b.v[2] - v[2] * b.v[1], v[2] * b.v[0] - v[0] * b.v[2],
          v[0] * b.v[1] - v[1] * b.v[0]};
}
inline double Vec3::sum() const {
  return v[0] + v[1] + v[2];
}
inline Vec3 Vec3::multiply(const Vec3& b) const {
  return Vec3(v[0] * b.v[0], v[1] * b.v[1], v[2] * b.v[2]);
}
inline Vec3& operator+=(Vec3& a, const Vec3& b) {
  a.v[0] += b.v[0];
  a.v[1] += b.v[1];
  a.v[2] += b.v[2];
  return a;
}
inline Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]);
}
inline Vec3& operator-=(Vec3& a, const Vec3& b) {
  a.v[0] -= b.v[0];
  a.v[1] -= b.v[1];
  a.v[2] -= b.v[2];
  return a;
}
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]);
}
inline Vec3& operator*=(Vec3& a, double b) {
  a.v[0] *= b;
  a.v[1] *= b;
  a.v[2] *= b;
  return a;
}
inline Vec3 operator*(const Vec3& a, const double b) {
  return Vec3(b * a.v[0], b * a.v[1], b * a.v[2]);
}
inline Vec3 operator*(const double b, const Vec3& a) {
  return Vec3(b * a.v[0], b * a.v[1], b * a.v[2]);
}
inline Vec3& operator/=(Vec3& a, double b) {
  a.v[0] /= b;
  a.v[1] /= b;
  a.v[2] /= b;
  return a;
}
inline Vec3 operator/(const Vec3& a, const double b) {
  return Vec3(a.v[0] / b, a.v[1] / b, a.v[2] / b);
}
inline double& Vec3::operator[](int index) { return v[index]; }
inline double Vec3::operator[](int index) const { return v[index]; }
inline double Vec3::norm() const { return std::sqrt(this->dot(*this)); }
// checks if a and b are on the same side of plane defined by normal
inline bool onSameSideOfPlane(const Vec3& a, const Vec3& b, const Vec3& normal) {
  return normal.dot(a) * normal.dot(b) > EPS;
}
// project a on b
inline Vec3 project(const Vec3& a, const Vec3& b) {
  return a.dot(b) / b.dot(b) * b;
}
// mirrors a over b
inline Vec3 mirrorOver(const Vec3& a, const Vec3& b) {
  return 2 * project(a, b) - a;
}

// A point on a triangle
struct TrianglePoint {
  const Triangle* triangle;
  Vec2 bary_coords;
};
class Voxel {
 public:
  Vec3 lo;
  Vec3 hi;
  Voxel() = default;
  Voxel(const Vec3& a, const Vec3& b);
  bool intersects(const Ray& r) const;
  // Clips triangle leaving side 'side' of plane
  void clip(const AxisPlane& plane, bool side);
  // Extends voxel to cover point p
  void cover(const Vec3& p);
  // Extends voxel to cover triangle t
  void cover(const Triangle& t);
  double area() const;
  // check that the point p is within the voxel (+EPS)
  bool isInside(const Vec3& p) const;
};
std::ostream& operator<<(std::ostream& out, const Voxel& a);
// Defines a point where ray intersects a triangle
struct RayIntersection {
  Vec2 bary_coords;
  double distance;
};
class Triangle {
 public:
  Vec3 p0, p1, p2;
  Triangle() = default;
  Triangle(const Vec3& p0, const Vec3& p1, const Vec3& p2);
  // Returns the point corresponding to barycentric coordinates coords
  Vec3 pointFromBary(const Vec2& coords) const;
  double area() const;
  void translate(const Vec3& v);
  // Finds non-parallel intersections between ray r and *this
  // distance == inf if no intersection was found
  // can return distance < 0 if the line r intersects
  // the triangle before the origin of r!
  RayIntersection getRayIntersection(const Ray& r) const;
};
std::ostream& operator<<(std::ostream& out, const Triangle& t);
class PlanePolygon {
 public:
  PlanePolygon(const Vec3& a, const Vec3& b, const Vec3& c);
  PlanePolygon(const Triangle& t);
  Voxel getBoundingBox() const;
  void intersect(const AxisPlane& plane, bool side);
  size_t size() const;
  Vec3& operator[](size_t index);
  Vec3 operator[](size_t index) const;

 private:
  std::vector<Vec3> points;
};
inline Vec3& PlanePolygon::operator[](size_t index) { return points[index]; }
inline Vec3 PlanePolygon::operator[](size_t index) const {
  return points[index];
}
// A class representing a clipped version of a triangle
//
// Clipped triangle is a shape S = t \cap clip_box
// where t is some triangle
// and clip_box is formed by the axis aligned splits
//
// Supports only queries about the bounding box of S
class ClipTriangle {
 public:
  ClipTriangle(const Triangle& triangle) : polygon(triangle) {
    box = polygon.getBoundingBox();
  }
  bool isAxisAligned(int axis) const;
  // returns subtrees the triangle should be added to
  std::pair<bool, bool> overlapsSides(const AxisPlane& plane, bool side) const;
  double max(int axis) const;
  double min(int axis) const;
  // Updates ClipTriangle so that side 'side' of plane is remains
  // The result should never be empty!
  void clip(const AxisPlane& plane, bool side);

 private:
  // Uses polygon representation of the clipped triangle to calculate
  // the bounding box.
  //
  // This is not necessary and is intended to be a temporary solution
  //
  // Notice that updating the bounding box is not trivial.
  // Consider following:
  // |\|->
  // | |
  // | |\   .
  // L-|-\  .
  // Here clip parallel to y axis induces changes the
  // minimal bounding box of the clip triangle both in
  // y and x directions
  PlanePolygon polygon;
  Voxel box;
};
inline bool ClipTriangle::isAxisAligned(int axis) const {
  return std::abs(box.hi[axis] - box.lo[axis]) < EPS;
}
class Ray {
 public:
  // p + sx;
  Vec3 origin;
  Vec3 direction;
  // inv_direction[i] = 1 / direction[i]
  // speeds up the triangle intersection calculations
  //
  // TODO
  // what if the user decides to modify direction? when to update this?
  Vec3 inv_direction;
  Ray(Vec3 origin, Vec3 direction);
};
std::ostream& operator<<(std::ostream& out, const Ray& a);
struct AxisPlane {
  int axis;
  double pos;
};
std::ostream& operator<<(std::ostream& out, const AxisPlane& a);
struct RayTriangleIntersection {
  size_t index;
  Vec2 bary_coords;
};
// index of the intersected triangle and barycentric coordinates
// returns {triangles.size(), Vec2(0)} if no intersection was found
RayTriangleIntersection firstRayTriangleIntersection(
    const std::vector<Triangle>& triangles, const Ray& r);
Voxel boundingBox(const std::vector<Triangle>& triangles);
inline bool pointOnSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
  return std::abs((a - p).norm() + (p - b).norm() - (a - b).norm()) < EPS;
}
// rotates a so that a[1] will point towards b
inline Vec3 rotateYTo(const Vec3& a, const Vec3& b) {
  Vec3 normal_1(-b[1], b[0], 0.0);
  normal_1 = normal_1 / normal_1.norm();
  Vec3 normal_2 = b.cross(normal_1);
  assert(std::abs(normal_2.norm() - 1) < EPS);
  return a[1] * b + a[0] * normal_1 + a[2] * normal_2;
}
template <typename Generator>
Vec3 uniformRandomHemispherePoint(Vec3 direction, Generator& g) {
  assert(direction.norm() > EPS);
  direction /= direction.norm();
  // Let X be the azimuth and Y the elevation of a uniformly distributed
  // random variable on the hemisphere.
  // Let U0 ~ U(0, 1)
  // Let U1 ~ U(0, 1)
  // P(y < a) = F_y(a) = (\int_{0}^{a} 2*pi * sin(x) \,dx)/2*pi = 1 - cos(a)
  // Now
  // Y ~ F_y^(-1)(U1) = acos(1 - U1)
  // --> Y ~ F_y^(-1)(1 - U1) = acos(U1)
  // X = 2 * PI * U0
  //
  // Let (p_x, p_y, p_z) be the point (X, Y) in cartesian coordinates
  // p_y = cos(Y) = U1
  // p_z = sin(X)*sin(Y) = sin(X) * sqrt(1 - U1^2)
  // p_x = cos(X)*sin(Y) = cos(X) * sqrt(1 - U1^2)

  std::uniform_real_distribution U(0.0, 1.0);
  double u_0 = 2 * kPi * U(g);
  double u_1 = U(g);
  Vec3 p(std::sin(u_0) * std::sqrt(1 - u_1 * u_1), u_1,
         std::cos(u_0) * std::sqrt(1 - u_1 * u_1));
  return rotateYTo(p, direction);
}
template <typename Generator>
Vec3 uniformRandomSpherePoint(Vec3 direction, Generator& g) {
  assert(direction.norm() > EPS);
  // TODO implement properly
  std::uniform_real_distribution dist(0.0, 1.0);
  if (dist(g) < 0.5) {
    return uniformRandomHemispherePoint(direction, g);
  }
  return uniformRandomHemispherePoint(-1.0 * direction, g);
}
// Samples points from a hemisphere with pdf
// (n + 1)/(2 * kPi) * direction.dot(v)^n
// see:
// https://graphics.cs.kuleuven.be/publications/Phong/.
template <typename Generator>
Vec3 cosineExponentRandomPoint(Vec3 direction, double exponent, Generator& g) {
  assert(direction.norm() > EPS);
  direction /= direction.norm();
  std::uniform_real_distribution U(0.0, 1.0);
  double e_1 = U(g);
  double e_2 = U(g);
  double cos_a = std::pow(e_1, 1.0 / (exponent + 1));
  double sin_a = std::sqrt(1.0 - cos_a*cos_a);
  Vec3 p(sin_a * std::cos(2.0 * kPi * e_2),
         cos_a,
         sin_a * std::sin(2.0 * kPi * e_2));
  return rotateYTo(p, direction);
}
// pdf of cosine exponent distribution facing 'direction'
// 0.0 outside the hemisphere!
inline double cosineExponentPdf(Vec3 v, Vec3 direction, double exponent) {
  v /= v.norm();
  direction.norm();
  if (v.dot(direction) < EPS) return 0.0;
  return (exponent + 1) / (2 * kPi) * std::pow(v.dot(direction), exponent);
}
#endif
