#ifndef RAYCASTER_GEOMETRY_H
#define RAYCASTER_GEOMETRY_H
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

const double INF = std::numeric_limits<double>::infinity();
const double EPS = 1e-10;

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
  friend Vec2 operator+(const Vec2& a, const Vec2& b);
  friend Vec2 operator-(const Vec2& a, const Vec2& b);
  friend Vec2 operator*(const Vec2& a, const double b);
  friend Vec2 operator*(const double b, const Vec2& a);
  double& operator[](int index);
  double operator[](int index) const;
  double norm() const;
  friend std::ostream& operator<<(std::ostream& out, const Vec2& a);

 private:
  double v[2];
};
inline double Vec2::dot(const Vec2& b) const { return v[0] * b.v[0] + v[1] * b.v[1]; }
inline Vec2 operator+(const Vec2& a, const Vec2& b) {
  return Vec2(a.v[0] + b.v[0], a.v[1] + b.v[1]);
}
inline Vec2 operator-(const Vec2& a, const Vec2& b) {
  return Vec2(a.v[0] - b.v[0], a.v[1] - b.v[1]);
}
inline Vec2 operator*(const Vec2& a, const double b) {
  return Vec2(b * a.v[0], b * a.v[1]);
}
inline Vec2 operator*(const double b, const Vec2& a) {
  return Vec2(b * a.v[0], b * a.v[1]);
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
  friend Vec3 operator+(const Vec3& a, const Vec3& b);
  friend Vec3 operator-(const Vec3& a, const Vec3& b);
  friend Vec3 operator*(const Vec3& a, const double b);
  friend Vec3 operator*(const double b, const Vec3& a);
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
inline Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]);
}
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]);
}
inline Vec3 operator*(const Vec3& a, const double b) {
  return Vec3(b * a.v[0], b * a.v[1], b * a.v[2]);
}
inline Vec3 operator*(const double b, const Vec3& a) {
  return Vec3(b * a.v[0], b * a.v[1], b * a.v[2]);
}
inline double& Vec3::operator[](int index) { return v[index]; }
inline double Vec3::operator[](int index) const { return v[index]; }
inline double Vec3::norm() const { return std::sqrt(this->dot(*this)); }
// A point on a triangle
struct TrianglePoint {
  Triangle* triangle;
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
  // Extends voxel to cover point t
  void cover(Triangle* t);
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
  PlanePolygon(const Triangle* t);
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
  Triangle* triangle;
  ClipTriangle(Triangle* triangle) : triangle(triangle), polygon(triangle) {
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
TrianglePoint firstRayTriangleIntersection(
    const std::vector<Triangle*>& triangles, const Ray& r);
Voxel boundingBox(const std::vector<Triangle*>& triangles);
#endif
