#include "material.h"

#include <cmath>
#include <tuple>

#include "lighting_utils.h"
// returns a random direction from importance sampling
Vec3 Material::importanceSample(const Vec3& normal, const Vec3& out_vector,
                                std::mt19937& mt) const {
  assert(std::abs(out_vector.norm() - 1.0) < EPS);
  assert(std::abs(normal.norm() - 1.0) < EPS);
  Vec3 direction;
  // probably could just return an error or any vector
  // but to keep the sampling separated from BRDF evaluation
  // try to return a reasonable direction
  if (diffuse.sum() + specular.sum() < EPS) {
    return sampleTransparentDiffuse(normal, mt);
  }
  double p_diffuse = diffuse.sum() / (diffuse.sum() + specular.sum());
  std::uniform_real_distribution uniform_random_0_1{0.0, 1.0};
  if (uniform_random_0_1(mt) < p_diffuse) {
    return sampleDiffuse(normal, mt);
  }
  return sampleSpecular(normal, out_vector, mt);
}
double Material::importanceSamplePdf(const Vec3& in_vector, const Vec3& normal,
                                     const Vec3& out_vector) const {
  assert(std::abs(out_vector.norm() - 1.0) < EPS);
  assert(std::abs(normal.norm() - 1.0) < EPS);
  double pdf = 0.0;
  if (diffuse.sum() + specular.sum() < EPS) {
    return transparentDiffusePdf(in_vector, normal);
  }
  double p_diffuse = diffuse.sum() / (diffuse.sum() + specular.sum());
  pdf += p_diffuse * diffusePdf(in_vector, normal);
  pdf += (1.0 - p_diffuse) * specularPdf(in_vector, normal, out_vector);
  return pdf;
}
// in_vector, normal and out_vector should be normalized
Vec3 Material::apply(const Vec3& in_vector, const Vec3& normal,
                     const Vec3& out_vector, const Vec3& in_color) const {
  assert(std::abs(out_vector.norm() - 1.0) < EPS);
  assert(std::abs(in_vector.norm() - 1.0) < EPS);
  assert(std::abs(normal.norm() - 1.0) < EPS);
  Vec3 total_brdf(0.0);
  total_brdf += diffuse.multiply(diffuseBrdf(in_vector, normal, out_vector));
  total_brdf += specular.multiply(specularBrdf(in_vector, normal, out_vector));
  return in_color.multiply(total_brdf) * std::abs(normal.dot(in_vector));
}
// samples from a cosine weighted distribution
Vec3 Material::sampleDiffuse(const Vec3& normal, std::mt19937& mt) const {
  double p_transparent = transparent.sum() / 3.0;
  std::uniform_real_distribution uniform_random_0_1{0.0, 1.0};
  if (uniform_random_0_1(mt) < p_transparent) {
    return sampleTransparentDiffuse(normal, mt);
  } else {
    return sampleOpaqueDiffuse(normal, mt);
  }
}
Vec3 Material::sampleTransparentDiffuse(const Vec3& normal,
                                        std::mt19937& mt) const {
  std::uniform_real_distribution uniform_random_0_1{0.0, 1.0};
  if (uniform_random_0_1(mt) < 0.5) {
    return cosineExponentRandomPoint(normal, 1.0, mt);
  } else {
    return cosineExponentRandomPoint(-1.0 * normal, 1.0, mt);
  }
}
Vec3 Material::sampleOpaqueDiffuse(const Vec3& normal, std::mt19937& mt) const {
  return cosineExponentRandomPoint(normal, 1.0, mt);
}
Vec3 Material::sampleSpecular(const Vec3& normal, const Vec3& out_vector,
                              std::mt19937& mt) const {
  double p_transparent = transparent.sum() / 3.0;
  std::uniform_real_distribution uniform_random_0_1{0.0, 1.0};
  if (uniform_random_0_1(mt) < p_transparent) {
    return sampleTransparentSpecular(normal, out_vector, mt);
  }
  return sampleOpaqueSpecular(normal, out_vector, mt);
}
Vec3 Material::sampleTransparentSpecular(const Vec3& normal,
                                         const Vec3& out_vector,
                                         std::mt19937& mt) const {
  Vec3 refraction;
  bool refraction_ok;
  std::tie(refraction, refraction_ok) =
      perfectRefraction(out_vector, normal, 1.0, index_of_refraction);
  Vec3 out_reflection = mirrorOver(out_vector, normal);
  // total internal reflection
  if (!refraction_ok) {
    return cosineExponentRandomPoint(out_reflection, specular_exp, mt);
  }
  double f = fresnelFactor(out_vector, normal, 1.0, index_of_refraction);
  std::uniform_real_distribution uniform_random_0_1{0.0, 1.0};
  if (uniform_random_0_1(mt) < f) {
    return cosineExponentRandomPoint(out_reflection, specular_exp, mt);
  }
  return cosineExponentRandomPoint(refraction, specular_exp, mt);
}
Vec3 Material::sampleOpaqueSpecular(const Vec3& normal, const Vec3& out_vector,
                                    std::mt19937& mt) const {
  Vec3 reflection = mirrorOver(out_vector, normal);
  return cosineExponentRandomPoint(reflection, specular_exp, mt);
}
double Material::diffusePdf(const Vec3& in_vector, const Vec3& normal) const {
  double p_transparent = transparent.sum() / 3.0;
  double pdf = 0.0;
  pdf += p_transparent * transparentDiffusePdf(in_vector, normal);
  pdf += (1.0 - p_transparent) * opaqueDiffusePdf(in_vector, normal);
  return pdf;
}
double Material::transparentDiffusePdf(const Vec3& in_vector,
                                       const Vec3& normal) const {
  double pdf = 0.0;
  pdf += 0.5 * cosineExponentPdf(in_vector, normal, 1.0);
  pdf += 0.5 * cosineExponentPdf(in_vector, -1.0 * normal, 1.0);
  return pdf;
}
double Material::opaqueDiffusePdf(const Vec3& in_vector,
                                  const Vec3& normal) const {
  return cosineExponentPdf(in_vector, normal, 1.0);
}
double Material::specularPdf(const Vec3& in_vector, const Vec3& normal,
                             const Vec3& out_vector) const {
  double p_transparent = transparent.sum() / 3.0;
  double pdf = 0.0;
  pdf += p_transparent * transparentSpecularPdf(in_vector, normal, out_vector);
  pdf +=
      (1.0 - p_transparent) * opaqueSpecularPdf(in_vector, normal, out_vector);
  return pdf;
}
double Material::transparentSpecularPdf(const Vec3& in_vector,
                                        const Vec3& normal,
                                        const Vec3& out_vector) const {
  Vec3 out_refraction;
  bool refraction_ok;
  std::tie(out_refraction, refraction_ok) =
      perfectRefraction(out_vector, normal, 1.0, index_of_refraction);
  Vec3 out_reflection = mirrorOver(out_vector, normal);
  // total internal reflection
  if (!refraction_ok) {
    return cosineExponentPdf(in_vector, out_reflection, specular_exp);
  }
  double pdf = 0.0;
  // note that here we have to use out_vector to estimate the incoming
  // ray fresnel factor because this is the pdf of the importance sampling
  // distribution which was generated with out_vector
  double f = fresnelFactor(out_vector, normal, 1.0, index_of_refraction);
  pdf += f * cosineExponentPdf(in_vector, out_reflection, specular_exp);
  pdf += (1.0 - f) * cosineExponentPdf(in_vector, out_refraction, specular_exp);
  return pdf;
}
double Material::opaqueSpecularPdf(const Vec3& in_vector, const Vec3& normal,
                                   const Vec3& out_vector) const {
  Vec3 reflection = mirrorOver(out_vector, normal);
  return cosineExponentPdf(in_vector, reflection, specular_exp);
}
Vec3 Material::diffuseBrdf(const Vec3& in_vector, const Vec3& normal,
                           const Vec3& out_vector) const {
  Vec3 out(0.0);
  out += transparent.multiply(transparentDiffuse());
  out += (Vec3(1.0) - transparent)
             .multiply(opaqueDiffuse(in_vector, normal, out_vector));
  return out;
}
Vec3 Material::transparentDiffuse() const { return 1.0 / (2.0 * kPi); }
Vec3 Material::opaqueDiffuse(const Vec3& in_vector, const Vec3& normal,
                             const Vec3& out_vector) const {
  if (onSameSideOfPlane(in_vector, out_vector, normal)) {
    return 1.0 / kPi;
  }
  return Vec3(0.0);
}
Vec3 Material::specularBrdf(const Vec3& in_vector, const Vec3& normal,
                            const Vec3& out_vector) const {
  Vec3 out(0.0);
  Vec3 tmp = transparentSpecular(in_vector, normal, out_vector);
  out += transparent.multiply(tmp);
  tmp = opaqueSpecular(in_vector, normal, out_vector);
  out += (Vec3(1.0) - transparent).multiply(tmp);
  return out;
}
Vec3 Material::transparentSpecular(const Vec3& in_vector, const Vec3& normal,
                                   const Vec3& out_vector) const {
  Vec3 in_refraction;
  bool refraction_ok;
  std::tie(in_refraction, refraction_ok) =
      perfectRefraction(in_vector, normal, 1.0, index_of_refraction);
  // total internal reflection
  if (!refraction_ok) {
    return phongSpecular(in_vector, normal, out_vector, specular_exp);
  }
  assert(std::abs(in_refraction.norm() - 1.0) < EPS);
  Vec3 out(0.0);
  double f = fresnelFactor(in_vector, normal, 1.0, index_of_refraction);
  out += f * phongSpecular(in_vector, normal, out_vector, specular_exp);
  assert(std::abs(in_vector.norm() - 1.0) < EPS);
  assert(std::abs(normal.norm() - 1.0) < EPS);
  // where the light would need to come to reflect like 'refraction'
  Vec3 imaginary_in = mirrorOver(in_refraction, normal);
  out +=
      (1.0 - f) * phongSpecular(imaginary_in, normal, out_vector, specular_exp);
  return out;
}
// Separate function so that adding additional grazing angle specular
// reflections would be possible in the future
Vec3 Material::opaqueSpecular(const Vec3& in_vector, const Vec3& normal,
                              const Vec3& out_vector) const {
  return phongSpecular(in_vector, normal, out_vector, specular_exp);
}
