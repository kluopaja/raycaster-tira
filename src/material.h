#ifndef RAYCASTER_MATERIAL_H
#define RAYCASTER_MATERIAL_H
#include <ostream>
#include <random>
#include <utility>

#include "geometry.h"

class Material {
 public:
  Vec3 diffuse = Vec3(0.0);
  Vec3 emitted = Vec3(0.0);
  Vec3 specular = Vec3(0.0);
  double specular_exp = 0.0;
  Vec3 transparent = Vec3(0.0);
  double index_of_refraction = 1.0;
  // Returns a direction vector from an imporatance sample distribution
  // This imporantance sample distribution is determined based only on the
  // material properties and out_vector but not, for example, light sources.
  //
  // `out_vector` and `normal` should be normalized
  Vec3 importanceSample(const Vec3& normal, const Vec3& out_vector,
                        std::mt19937& mt) const;
  // Pdf of the in_vector sampled from importanceSample distribution.
  // Used with the Monte Carlo integration.
  //
  // `in_vector`, `normal` and `out_vector` should be normalized.
  double importanceSamplePdf(const Vec3& in_vector, const Vec3& normal,
                             const Vec3& out_vector) const;
  // Applies material to light of color `in_color` coming from `in_vector`
  // and exiting to `out_vector`
  //
  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 apply(const Vec3& in_vector, const Vec3& normal, const Vec3& out_vector,
             const Vec3& in_color) const;

 private:
  // For importance sampling
  // Sampling functions for different types of lighting
  //
  // `normal` should be normalized
  Vec3 sampleDiffuse(const Vec3& normal, std::mt19937& mt) const;
  // `normal` should be normalized
  Vec3 sampleTransparentDiffuse(const Vec3& normal, std::mt19937& mt) const;
  // `normal` should be normalized
  Vec3 sampleOpaqueDiffuse(const Vec3& normal, std::mt19937& mt) const;
  // `normal` and `out_vector` should be normalized
  Vec3 sampleSpecular(const Vec3& normal, const Vec3& out_vector,
                      std::mt19937& mt) const;
  // `normal` and `out_vector` should be normalized
  Vec3 sampleTransparentSpecular(const Vec3& normal, const Vec3& out_vector,
                                 std::mt19937& mt) const;
  // `normal` and `out_vector` should be normalized
  Vec3 sampleOpaqueSpecular(const Vec3& normal, const Vec3& out_vector,
                            std::mt19937& mt) const;

  // For pdf calculations
  //
  // Pdfs for different types of lighting
  //
  // `in_vector` and `normal` should be normalized
  double diffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  // `in_vector` and `normal` should be normalized
  double transparentDiffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  // `in_vector` and `normal` should be normalized
  double opaqueDiffusePdf(const Vec3& in_vector, const Vec3& normal) const;
  // `in_vector`,  `normal` and `out_vector` should be normalized
  double specularPdf(const Vec3& in_vector, const Vec3& normal,
                     const Vec3& out_vector) const;
  // `in_vector`,  `normal` and `out_vector` should be normalized
  double transparentSpecularPdf(const Vec3& in_vector, const Vec3& normal,
                                const Vec3& out_vector) const;
  // `in_vector`,  `normal` and `out_vector` should be normalized
  double opaqueSpecularPdf(const Vec3& in_vector, const Vec3& normal,
                           const Vec3& out_vector) const;
  // For applying the material
  //
  // Applying different types of lighting
  //
  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 diffuseBrdf(const Vec3& in_vector, const Vec3& normal,
                   const Vec3& out_vector) const;
  Vec3 transparentDiffuse() const;
  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 opaqueDiffuse(const Vec3& in_vector, const Vec3& normal,
                     const Vec3& out_vector) const;

  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 specularBrdf(const Vec3& in_vector, const Vec3& normal,
                    const Vec3& out_vector) const;
  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 transparentSpecular(const Vec3& in_vector, const Vec3& normal,
                           const Vec3& out_vector) const;
  // `in_vector`, `normal` and `out_vector` should be normalized
  Vec3 opaqueSpecular(const Vec3& in_vector, const Vec3& normal,
                      const Vec3& out_vector) const;
};
std::ostream& operator<<(std::ostream& out, const Material& m);
#endif
