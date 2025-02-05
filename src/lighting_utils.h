#ifndef RAYCASTER_LIGHTING_UTILS_H
#define RAYCASTER_LIGHTING_UTILS_H
#include "geometry.h"
// calculates the specular component of Phong lighting
// returns Vec3(0.0) if in_vector and out_vector are on
// different sides of plane defined by 'normal'
// all input vectors should be normalized
inline Vec3 phongSpecular(const Vec3& in_vector, const Vec3& normal,
                          const Vec3& out_vector, double exponent) {
  assert(std::abs(in_vector.norm() - 1) < EPS);
  assert(std::abs(normal.norm() - 1) < EPS);
  assert(std::abs(out_vector.norm() - 1) < EPS);
  if (!onSameSideOfPlane(in_vector, out_vector, normal)) {
    return Vec3(0.0);
  }
  // To follow the conservation of energy we need to
  // multiply by this
  double normalization = (exponent + 1.0) / (2 * kPi);
  double cos_alpha =
      std::max(out_vector.dot(mirrorOver(in_vector, normal)), 0.0);
  double cos_theta = std::abs(normal.dot(out_vector));
  return normalization * std::pow(cos_alpha, exponent) / cos_theta;
}
// eta_1 is the refractive index of material on the side of normal
// eta_2 is the refractive index of material on the other side
// returns {Vec3(0.0), 0} if there was a total internal reflection
// otherwise returns {v, 1} where v is the refracted ray
//
// Assumes that `in_vector` and `normal` are normalized
inline std::pair<Vec3, bool> perfectRefraction(Vec3 in_vector, Vec3 normal,
                                               double eta_1, double eta_2) {
  assert(std::abs(in_vector.norm() - 1.0) < EPS);
  assert(std::abs(normal.norm() - 1.0) < EPS);
  // internally the function assumes that the eta_1 is on the side of
  // plane where in_vector is
  // also that normal is on the same side as in_vector
  if (normal.dot(in_vector) < -EPS) {
    normal = -1.0 * normal;
    std::swap(eta_1, eta_2);
  }
  // Assume n and in_vector are normalized
  // Let w1 = proj_n(in_vector)
  // and w2 = proj_n(out_vector)
  // Let u1 = in_vector - w1
  // and u2 = out_vector - w2
  // Let t1 be the angle between in_vector and n
  // and t2 be the angle between out_vector and -n
  //
  // sin(t1) / sin(t2) = eta_2 / eta_1
  // --> sin(t2) = sin(t1) * eta_1 / eta_2
  // --> |u2| = |u1| * eta_1 / eta_2
  // --> u2 = -u1 * eta_1 / eta_2
  //
  // cos(t2) = sqrt(1 - sin(t2)^2)
  // |w2| = cos(t2) = sqrt(1 - sin(t2)^2)
  // w2 = -n * sqrt(1 - dot(u2, u2))
  Vec3 w1 = project(in_vector, normal);
  Vec3 u1 = in_vector - w1;
  Vec3 u2 = -1.0 * u1 * eta_1 / eta_2;
  double u2_dot = u2.dot(u2);
  if (u2_dot > 1.0 - EPS) {
    return {Vec3(0.0), 0};
  }
  Vec3 w2 = -1.0 * normal * std::sqrt(1.0 - u2_dot);
  return {(u2 + w2), 1};
}
// calculates fresnel factor using the Schlick's approximation
// eta_1 is the refractive index of the side on the side of normal
// eta_2 is the refractive index of the other side
// In Schlick's approximation these are identical!
// All input vectors should be normalized
inline double fresnelFactor(const Vec3& in_vector, const Vec3& normal,
                            double eta_1, double eta_2) {
  // note that swapping \eta_1 and \eta_2 would just swap the
  // sign of tmp. This would not be present in tmp * tmp
  double tmp = (eta_1 - eta_2) / (eta_1 + eta_2);
  double r_0 = tmp * tmp;
  return r_0 +
         (1.0 - r_0) * std::pow(1.0 - std::abs(normal.dot(in_vector)), 5.0);
}
#endif
