#ifndef RAYCASTER_RAYCASTER_H
#define RAYCASTER_RAYCASTER_H
#include <random>
#include <vector>

#include "geometry.h"
#include "image.h"
#include "kd_tree.h"
#include "material.h"
enum NormalType {
  kSmooth,
  kRough,
};
// sampling scheme for shooting rays
// kImportanceSampling is the default
// kUniformSphere might be useful for debugging
enum SamplingScheme {
  kImportanceSampling,
  kUniformSphere,
};
struct SceneTriangle {
  Triangle triangle;
  Vec3 normals[3];
  Material material;
  // returns a normalized interpolated triangle normal for the triangle point
  // bary_coords
  Vec3 normalAt(const Vec2& bary_coords);
};
class Model {
 public:
  Model() = default;
  Model(Model&& a) = default;
  Model& operator=(Model&& a) = default;
  void concatenate(const Model& a);
  void translate(const Vec3& v);
  Tree buildKdTree(double k_t, double k_i);
  Vector<SceneTriangle> scene_triangles;
  Vector<Material> materials;
};
class PointLight {
 public:
  PointLight(const Vec3& position, const Vec3& color);
  Vec3 position;
  Vec3 color;
};
class Camera {
 public:
  // 0 < x_fov, y_fov < PI
  // position + front_vec will be the middle point of the
  // image plane
  // only things in the front of the image plane will be rendered
  Camera(const Vec3& position, const Vec3& front_vec, const Vec3& up_vec,
         double x_fov, double y_fov);

  // returns a ray 'r' starting from image plane coordinates (x, y)
  // and extending outside the image plane.
  // A line parallel to 'r' will intersect 'position'
  Ray rayFromImagePlane(double x, double y);

 private:
  // position of the camera
  Vec3 position;
  // image_position = (left top of the image plane)
  Vec3 image_position;
  // image_right_vec
  // = (right top of the image plane) - (left top of the image plane)
  Vec3 image_right_vec;
  // image_down_vec
  // = (left bottom of the image plane) - (left top of the image plane)
  Vec3 image_down_vec;
};
// Class representing the environment light, that is, that
// a ray receives if it doesn't intersect with any triangle
class EnvironmentLight {
  public:
   // All directections emit light with intensity `this->color`
   void setUniform();
   // L from direction v will have the
   // radiance of `this->color` * pow(max(0, cos(alpha)), `exponent`)
   // where alpha is the angle between `direction` and v
   void setDirected(const Vec3& direction, double exponent);
   void setColor(const Vec3& color);
   // Calculates the radiance of light coming from direciton `d`
   Vec3 colorAtDirection(Vec3 direction);
  private:
   // speeds up the queries a bit if this is handled separately
   bool is_directed = 0;
   Vec3 color = Vec3(0.0);
   Vec3 direction = Vec3(0.0);
   // should be >= 0
   double cosine_exp = 0;
};
// std::pow(base, 0) always returns 1
inline Vec3 EnvironmentLight::colorAtDirection(Vec3 direction) {
  if(!is_directed) {
    return color;
  }
  else {
    direction /= direction.norm();
    double cos_a = std::max(0.0, direction.dot(this->direction));
    return color.multiply(std::pow(cos_a, cosine_exp));
  }
}
// For parallel rendering
struct ThreadResources {
  std::mt19937 thread_mt;
  Vector<std::pair<int, int> > pixels;
};
class Scene {
 public:
  Scene(const Camera& c);
  bool addModelFromFile(const std::string& file, const Vec3& position,
                        NormalType normal_type);
  void addPointLight(const Vec3& position, const Vec3& color);
  // Uniform environment light from every direction
  void setEnvironmentLightUniform();
  // Environment light from direction v will have the
  // radiance of environment_light_color * cos(alpha)^exponent
  // where environment_light_color can be set with
  // setEnvironmentLightColor and alpha is the angle between
  // `direction` and v
  void setEnvironmentLightDirected(const Vec3& direction, double exponent);
  void setEnvironmentLightColor(const Vec3& color);
  void setSamplingScheme(SamplingScheme s);
  Image render(int x_resolution, int y_resolution, int n_rays_per_pixel,
               int max_recursion_depth, int n_recursion_rays);

 private:
  Vec3 renderPixel(int x_pixel, int y_pixel, int x_resolution, int y_resolution,
                   int n_rays_per_pixel, std::mt19937& thread_mt_19937);
  Vec3 renderImagePlanePoint(double x, double y, std::mt19937& thread_mt_19937);
  Vec3 castRay(const Ray& r, int recursion_depth,
               std::mt19937& thread_mt_19937);
  Vec3 totalPointLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                            const Material& material);
  Vec3 pointLightColor(const Vec3& point, const Vec3& normal,
                       const Vec3& out_vector, const Material& material,
                       const PointLight& point_light);
  Vec3 indirectLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                          const Material& material, int depth,
                          std::mt19937& thread_mt_19937);
  Vec3 sampleIndirectLight(const Vec3& point, const Vec3& normal,
                           const Vec3& out_vector, const Material& material,
                           int recursion_depth, std::mt19937& thread_mt_19937);
  Vec3 environmentLightColor(const Ray& r);
  EnvironmentLight environment_light;
  Model model;
  Vector<PointLight> point_lights;
  Camera camera;
  Tree kd_tree;
  std::mt19937 mt_19937;
  std::uniform_real_distribution<double> subpixel_sample_distribution;
  int max_recursion_depth;
  int n_recursion_rays;
  SamplingScheme sampling_scheme;
};
#endif
