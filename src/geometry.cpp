#include "geometry.h"
Vec2::Vec2(double x) : v{x, x} {}
Vec2::Vec2(double x, double y) : v{x, y} {}
std::ostream& operator<<(std::ostream& out, const Vec2& a) {
  out << "Vec2(" << a[0] << ", " << a[1] << ")";
  return out;
}
Vec3::Vec3(double x) : v{x, x, x} {}
Vec3::Vec3(double x, double y, double z) : v{x, y, z} {}
std::ostream& operator<<(std::ostream& out, const Vec3& a) {
  out << "Vec3(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
  return out;
}
std::istream& operator>>(std::istream& in, Vec3& a) {
  in >> a[0] >> a[1] >> a[2];
  return in;
}

Voxel::Voxel(const Vec3& lo, const Vec3& hi) : lo(lo), hi(hi) {}
// Implements a fast voxel ray intersection with the correct handling of
// special cases
//
// ( from here https://dl.acm.org/doi/10.1145/1198555.1198748 )
//
// ( Note that this algorithm doesn't always have consistent behaviour
// at the edges of the voxel! )
//
// A ray r intersect a voxel v iff there is a parameter t > 0
// s.t. r.origin + t * r.direction is inside the voxel v.
//
// A point p is inside a voxel v exactly when:
//      v.lo[0] < p[0] < v.hi[0]
// and  v.lo[1] < p[1] < v.hi[1]
// and  v.lo[2] < p[2] < v.hi[2]
//
// So we want to find a value t such that.
//      v.lo[0] < r.origin[0] + t * r.direction[0]  < v.hi[0]
//      v.lo[1] < r.origin[1] + t * r.direction[1]  < v.hi[1]
//      v.lo[2] < r.origin[2] + t * r.direction[2]  < v.hi[2]
//
// Any one of these can be solved simply:
//
// v.lo[i] < r.origin[i] + t * r.direction[i] < v.hi[i]
// (v.lo[i] - r.origin[i]) / r.direction[i] < t
//     < (v.hi[i] - r.origin[i]) / r.direction[i]
//
// (inversing the inequalities if r.direction[i] < 0)
//
// The algorithm handles these intervals one by one
// always maintaining the intersection [t_min, t_max] of
// the so far met intervals. Whenever possible,
// the algorithm checks whether the interval is empty
//
// While in principle the process is very simple there
// are some special cases when r.direction[i] == 0.0
// or r.direction[i] == -0.0
//
// While it holds that 0.0 == -0.0, 0.0 and -0.0
// still behave differently since 1/-0.0 == -inf < inf == 1/0.0.
// Therefore it is important that we compare divx >= 0 and
// not r.direction[i] >= 0 !
//
// If r.direction[i] == +-0.0, the [t_min, t_max]
// should be updated based on whether
// lo[i] < r.origin[i] < hi[i] holds
// if this inequality does not hold, the algorithm
// will set either (t_min = inf; t_max = inf)
// or (t_min = -inf; t_max = -inf)
// both of which define an interval that evaluates
// to 0 at the end of the algorithm!
//
// Other special case is r.direction[i] == +-0.0
// and (lo[i] == r.origin[i] or hi[i] == r.origin[i])
//
// Now we get for example:
//     t_min = (lo[0] - r.origin[0]) * divx
//           = 0.0*inf = NaN
//
// In these cases, the algorithm doesn't handle the
// dimensions consistently:
//     If this occurs for dimension 0, the algorithm returns false
//         (NaN comparisons will evaluate to false)
//     In later steps, these NaN values mean that the
//     corresponding limit (ty_min, ty_max, tz_min, tz_max)
//     will not limit t_min or t_max. So effectively in these cases
//     the intersection with the voxel border will be treated
//     as intersection with the voxel!

bool Voxel::intersects(const Ray& r) const {
  double t_min, t_max;
  double ty_min, ty_max;
  if (r.getInvDirection()[0] >= 0) {
    t_min = (lo[0] - r.getOrigin()[0]) * r.getInvDirection()[0];
    t_max = (hi[0] - r.getOrigin()[0]) * r.getInvDirection()[0];
  } else {
    t_min = (hi[0] - r.getOrigin()[0]) * r.getInvDirection()[0];
    t_max = (lo[0] - r.getOrigin()[0]) * r.getInvDirection()[0];
  }
  if (r.getInvDirection()[1] >= 0) {
    ty_min = (lo[1] - r.getOrigin()[1]) * r.getInvDirection()[1];
    ty_max = (hi[1] - r.getOrigin()[1]) * r.getInvDirection()[1];
  } else {
    ty_min = (hi[1] - r.getOrigin()[1]) * r.getInvDirection()[1];
    ty_max = (lo[1] - r.getOrigin()[1]) * r.getInvDirection()[1];
  }
  // check if the new intervals doesn't overlap with the first one
  if ((t_min > ty_max) || (ty_min > t_max)) {
    return false;
  }
  // if clauses and not std::min are used here because this way
  // it is easier to see happens with NaNs
  if (ty_min > t_min) {
    t_min = ty_min;
  }
  if (ty_max < t_max) {
    t_max = ty_max;
  }
  double tz_min, tz_max;
  if (r.getInvDirection()[2] >= 0) {
    tz_min = (lo[2] - r.getOrigin()[2]) * r.getInvDirection()[2];
    tz_max = (hi[2] - r.getOrigin()[2]) * r.getInvDirection()[2];
  } else {
    tz_min = (hi[2] - r.getOrigin()[2]) * r.getInvDirection()[2];
    tz_max = (lo[2] - r.getOrigin()[2]) * r.getInvDirection()[2];
  }
  if ((t_min > tz_max) || (tz_min > t_max)) {
    return false;
  }
  if (tz_min > t_min) {
    t_min = tz_min;
  }
  if (tz_max < t_max) {
    t_max = tz_max;
  }
  // We need to remember to handle the case t_min == inf
  return ((t_min < INF) && (t_max > 0));
}
void Voxel::clip(const AxisPlane& plane, bool side) {
  // shouldn't clip away the whole voxel!
  assert(lo[plane.axis] - EPS < plane.pos);
  assert(plane.pos <= hi[plane.axis] + EPS);
  if (side == 0) {
    hi[plane.axis] = plane.pos;
  } else {
    lo[plane.axis] = plane.pos;
  }
}
// extends Voxel to cover point p
void Voxel::cover(const Vec3& p) {
  for (int i = 0; i < 3; ++i) {
    lo[i] = std::min(lo[i], p[i]);
    hi[i] = std::max(hi[i], p[i]);
  }
}
// extends Voxel to cover triangle t
void Voxel::cover(const Triangle& t) {
  cover(t.p0);
  cover(t.p1);
  cover(t.p2);
};
double Voxel::area() const {
  return 2.0 * ((hi[0] - lo[0]) * (hi[1] - lo[1]) +
                (hi[1] - lo[1]) * (hi[2] - lo[2]) +
                (hi[2] - lo[2]) * (hi[0] - lo[0]));
}
bool Voxel::isInside(const Vec3& p) const {
  return lo[0] - EPS < p[0] && p[0] < hi[0] + EPS && lo[1] - EPS < p[1] &&
         p[1] < hi[1] + EPS && lo[2] - EPS < p[2] && p[2] < hi[2] + EPS;
}

std::ostream& operator<<(std::ostream& out, const Voxel& a) {
  out << "Voxel(\n" << a.lo << ",\n" << a.hi << ")";
  return out;
}
Triangle::Triangle(const Vec3& p0, const Vec3& p1, const Vec3& p2)
    : p0(p0), p1(p1), p2(p2) {}
Vec3 Triangle::pointFromBary(const Vec2& coords) const {
  assert(coords[0] + coords[1] < 1 + EPS);
  return (1.0 - coords[0] - coords[1]) * p0 + coords[0] * p1 + coords[1] * p2;
}
double Triangle::area() const { return (p1 - p0).cross(p2 - p0).norm() / 2.0; }
void Triangle::translate(const Vec3& v) {
  p0 = p0 + v;
  p1 = p1 + v;
  p2 = p2 + v;
}
// For explanation see
// https://cadxfem.org/inf/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
// also detects intersections to the side of the triangle
// and EPS outside the triangle
inline RayIntersection Triangle::getRayIntersection(const Ray& r) const {
  Vec3 edge1 = p1 - p0;
  Vec3 edge2 = p2 - p0;
  Vec3 pvec = r.getDirection().cross(edge2);
  double det = edge1.dot(pvec);
  // ray parallel to the triangle
  if (det > -EPS && det < EPS) {
    return {{}, std::numeric_limits<double>::infinity()};
  }
  double inv_det = 1.0 / det;
  Vec3 tvec = r.getOrigin() - p0;
  RayIntersection result;
  result.bary_coords[0] = tvec.dot(pvec) * inv_det;
  if (result.bary_coords[0] < -EPS || result.bary_coords[0] > 1.0 + EPS) {
    return {{}, std::numeric_limits<double>::infinity()};
  }
  Vec3 qvec = tvec.cross(edge1);
  result.bary_coords[1] = r.getDirection().dot(qvec) * inv_det;
  if (result.bary_coords[1] < -EPS ||
      result.bary_coords[0] + result.bary_coords[1] > 1.0 + EPS) {
    return {{}, std::numeric_limits<double>::infinity()};
  }
  result.distance = edge2.dot(qvec) * inv_det;
  return result;
}
std::ostream& operator<<(std::ostream& out, const Triangle& t) {
  out << "Triangle(\n" << t.p0 << ",\n" << t.p1 << ",\n" << t.p2 << ")";
  return out;
};
// A class representing a convex planar polygon in 3d space
PlanePolygon::PlanePolygon(const Vec3& a, const Vec3& b, const Vec3& c)
    : points{a, b, c} {
  assert((c - a).cross(b - a).norm() > EPS && "Degenerate triangle");
}
PlanePolygon::PlanePolygon(const Triangle& t)
    : PlanePolygon(t.p0, t.p1, t.p2) {}
Voxel PlanePolygon::getBoundingBox() const {
  Voxel bounds(INF, -INF);
  for (size_t i = 0; i < points.size(); ++i) {
    bounds.cover(points[i]);
  }
  return bounds;
}
void PlanePolygon::intersect(const AxisPlane& plane, bool side) {
  Vector<Vec3> new_points;
  for (size_t i = 0; i < points.size(); ++i) {
    bool retain = 0;
    if (side == 0) {
      if (points[i][plane.axis] < plane.pos + EPS) {
        retain = 1;
      }
    } else {
      if (points[i][plane.axis] > plane.pos - EPS) {
        retain = 1;
      }
    }
    if (retain) {
      new_points.pushBack(points[i]);
    }
    size_t nxt_i = i + 1;
    if (nxt_i == points.size()) nxt_i = 0;

    Vec3 a = points[i];
    Vec3 b = points[nxt_i];
    double side_length = b[plane.axis] - a[plane.axis];
    if (side_length < 0) {
      std::swap(a, b);
      side_length = -side_length;
    }
    // check if the plane intersect the segment
    if (a[plane.axis] + EPS < plane.pos && plane.pos < b[plane.axis] - EPS) {
      double t = (plane.pos - a[plane.axis]) / side_length;
      new_points.pushBack((1.0 - t) * a + t * b);
    }
  }
  points = new_points;
}
size_t PlanePolygon::size() const { return points.size(); }
// Returns whether triangle has positive area overlap
// with side 0 and 1 of plane
// 'side' is the side which the plane belongs to
// (side which a triangle lying on 'plane' should be added to)
std::pair<bool, bool> ClipTriangle::overlapsSides(const AxisPlane& plane,
                                                  bool side) const {
  if (isAxisAligned(plane.axis)) {
    if (box.hi[plane.axis] < plane.pos - EPS) {
      return {1, 0};
    }
    if (box.lo[plane.axis] > plane.pos + EPS) {
      return {0, 1};
    }
    if (side == 0) {
      return {1, 0};
    }
    return {0, 1};
  }
  if (box.hi[plane.axis] < plane.pos + EPS) {
    return {1, 0};
  }
  if (box.lo[plane.axis] > plane.pos - EPS) {
    return {0, 1};
  }
  return {1, 1};
}
double ClipTriangle::max(int axis) const { return box.hi[axis]; }
double ClipTriangle::min(int axis) const { return box.lo[axis]; }
void ClipTriangle::clip(const AxisPlane& plane, bool side) {
  assert(side != 0 || (plane.pos > min(plane.axis) - EPS));
  assert(side != 1 || (plane.pos < max(plane.axis) + EPS));
  polygon.intersect(plane, side);
  assert(polygon.size() > 1);
  box = polygon.getBoundingBox();
}
Ray::Ray(Vec3 origin, Vec3 direction) : origin(origin), direction(direction) {
  assert(direction.norm() > EPS);
  for (int i = 0; i < 3; ++i) {
    inv_direction[i] = 1.0 / direction[i];
  }
}
std::ostream& operator<<(std::ostream& out, const Ray& a) {
  out << "Ray(\n" << a.getOrigin() << ",\n" << a.getDirection() << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const AxisPlane& a) {
  out << "AxisPlane(\n" << a.axis << ",\n" << a.pos << ")";
  return out;
}
RayTriangleIntersection firstRayTriangleIntersection(
    const Vector<Triangle>& triangles, const Ray& r) {
  RayTriangleIntersection closest_point = {triangles.size(), {}};
  double closest_distance = INF;
  for (size_t i = 0; i < triangles.size(); ++i) {
    RayIntersection intersection = triangles[i].getRayIntersection(r);
    if (intersection.distance > 0 && intersection.distance < closest_distance) {
      closest_distance = intersection.distance;
      closest_point = {i, intersection.bary_coords};
    }
  }
  return closest_point;
}
Voxel boundingBox(const Vector<Triangle>& triangles) {
  if (triangles.size() == 0) {
    return {Vec3(0), Vec3(0)};
  }
  Voxel box = {Vec3(std::numeric_limits<double>::infinity()),
               Vec3(-std::numeric_limits<double>::infinity())};
  for (size_t i = 0; i < triangles.size(); ++i) {
    box.cover(triangles[i]);
  }
  return box;
}
