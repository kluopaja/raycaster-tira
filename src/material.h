#ifndef RAYCASTER_MATERIAL_H
#define RAYCASTER_MATERIAL_H
#include <random>
#include <utility>
#include <ostream>

#include "geometry.h"

class Material {
 public:
  Vec3 diffuse = Vec3(0.0);
  Vec3 emitted = Vec3(0.0);
  Vec3 specular = Vec3(0.0);
  double specular_exp = 0.0;
  Vec3 transparent= Vec3(0.0);
  double index_of_refraction = 1.0;
  // returns a direction vector from an imporatance sample distribution
  // This imporantance sample distribution is determined based only on the
  // material properties and out_vector but not, for example, light sources
  // out_vector and normal should be normalized
  Vec3 importanceSample(const Vec3& normal, const Vec3& out_vector,
                        std::mt19937& mt) const;
  // pdf of the in_vector in importance sample distribution
  // used with the Monte Carlo integration
  double importanceSamplePdf(const Vec3& in_vector, const Vec3& normal,
                             const Vec3& out_vector) const;
  // applies material on in_color
  // in_vector, normal and out_vector should be normalized
  Vec3 apply(const Vec3& in_vector, const Vec3& normal, const Vec3& out_vector,
             const Vec3& in_color) const;

 private:
  // for importance sampling
  Vec3 sampleDiffuse(const Vec3& normal, std::mt19937& mt) const;
  Vec3 sampleTransparentDiffuse(const Vec3& normal, std::mt19937& mt) const;
  Vec3 sampleOpaqueDiffuse(const Vec3& normal, std::mt19937& mt) const;
  Vec3 sampleSpecular(const Vec3& normal, const Vec3& out_vector,
                      std::mt19937& mt) const;
  Vec3 sampleTransparentSpecular(const Vec3& normal, const Vec3& out_vector,
                                 std::mt19937& mt) const;
  Vec3 sampleOpaqueSpecular(const Vec3& normal, const Vec3& out_vector,
                            std::mt19937& mt) const;
  // for pdf calculations
  double diffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  double transparentDiffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  double opaqueDiffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  double specularPdf(const Vec3& in_vector, const Vec3& normal,
                     const Vec3& out_vector) const;
  double transparentSpecularPdf(const Vec3& in_vector, const Vec3& normal,
                                const Vec3& out_vector) const;
  double opaqueSpecularPdf(const Vec3& in_vector, const Vec3& normal,
                           const Vec3& out_vector) const;
  // for applying the material
  Vec3 diffuseBrdf(const Vec3& in_vector, const Vec3& normal,
                   const Vec3& out_vector) const;
  Vec3 transparentDiffuse() const;
  Vec3 opaqueDiffuse(const Vec3& in_vector, const Vec3& normal,
                     const Vec3& out_vector) const;

  Vec3 specularBrdf(const Vec3& in_vector, const Vec3& normal,
                    const Vec3& out_vector) const;
  Vec3 transparentSpecular(const Vec3& in_vector, const Vec3& normal,
                           const Vec3& out_vector) const;
  Vec3 opaqueSpecular(const Vec3& in_vector, const Vec3& normal,
                      const Vec3& out_vector) const;
};
std::ostream& operator<<(std::ostream& out, const Material& m);
#endif
