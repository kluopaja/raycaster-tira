#ifndef RAYCASTER_RAYCASTER_H
#define RAYCASTER_RAYCASTER_H
#include <random>

#include "geometry.h"
#include "image.h"
#include "kd_tree.h"
#include "material.h"
enum class NormalType {
  kSmooth,
  kRough,
  N_ENUM_VALUES,
};
enum class EnvironmentLightType {
  kUniform,
  kDirected,
  N_ENUM_VALUES,
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
// Class for storing the 3D models for rendering
class Model {
 public:
  Model() = default;
  Model(Model&& a) = default;
  Model& operator=(Model&& a) = default;
  // Adds the content of `a` to `*this`
  void concatenate(const Model& a);
  // Translates the models coorinates by `v`
  void translate(const Vec3& v);
  // Returns a kd tree containing the triangles of the object
  //
  // To save memory (a single triangle can be stored in many nodes in
  // the kd tree), the kd tree stores pointers to
  // this->scene_triangles so scene_triangles should not
  // be modified after building the kd tree
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
// Class representing the viewpoint of the rendered image
class Camera {
 public:
  // 0 < x_fov, y_fov < PI
  // position + front_vec will be the middle point of the
  // image plane
  // only things in the front of the image plane will be rendered
  Camera(const Vec3& position, const Vec3& front_vec, const Vec3& up_vec,
         double x_fov, double y_fov);

  // Returns a ray 'r' starting from image plane coordinates (x, y)
  // and extending outside the image plane.
  // A line parallel to 'r' will intersect 'position'
  //
  // Assumes 0 <= x, y <= 1
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
  // The random number generator used by the thread
  std::mt19937 thread_mt;
  // The pixels that the thread will render
  Vector<std::pair<int, int> > pixels;
};
// Class storing a representation of a rendering scene
// Can be used to render a scene to an Image object
class Scene {
 public:
  Scene(const Camera& c);
  // Imports a model from file `file`, translates it to `position`
  // and adds it to the current scene
  //
  // If the model includes definitions for normals, then `normal_type`
  // is ignored. Otherwise `normal_type` determines the type
  // of normals generated by the importer.
  bool addModelFromFile(const std::string& file, const Vec3& position,
                        NormalType normal_type);
  // Adds a point light of `color` to `position`.
  // `color` should not contain negative elements.
  void addPointLight(const Vec3& position, const Vec3& color);
  // Uniform environment light from every direction
  void setEnvironmentLightUniform();
  // Environment light from direction v will have the
  // radiance of environment_light_color * cos(alpha)^exponent
  // where environment_light_color can be set with
  // setEnvironmentLightColor and alpha is the angle between
  // `direction` and v
  void setEnvironmentLightDirected(const Vec3& direction, double exponent);
  // `color` should not contain negative elements
  void setEnvironmentLightColor(const Vec3& color);
  // Determines the sampling method used to shoot recursive
  // rays after the first ray intersection.
  void setSamplingScheme(SamplingScheme s);
  // Sets cost estimates used in kd tree building
  // `new_k_t` and `new_k_i` should be positive
  void setCosts(double new_k_t, double new_k_i);
  // Renders the scene.
  //
  // `n_rays_per_pixel` is the number of rays sampled from every pixel
  //    if `n_rays_per_pixel` == 1, then the middle point is sampled.
  //    Otherwise `n_rays_per_pixel` rays will be sampled uniformly random
  //
  // `max_recursion_depth` is the maximum number of rays on the path
  // from the camera to any sampled intersection point in the scene.
  // Note that doesn't include the possible point light sampling rays.
  //
  // `n_recursion_rays` is the number of rays shot from an ray intersection
  // point excluding rays to point lights.
  //
  // Assumes all parameters are positive
  Image render(int x_resolution, int y_resolution, int n_rays_per_pixel,
               int max_recursion_depth, int n_recursion_rays);

 private:
  // Returns an estimate for the color of pixel (`x_pixel`, `y_pixel`)
  // `thread_mt_19937` should be accessed only within one thread
  Vec3 renderPixel(int x_pixel, int y_pixel, int x_resolution, int y_resolution,
                   int n_rays_per_pixel, std::mt19937& thread_mt_19937);
  // Returns an estimate for the image plane point (x, y)
  // Assumes 0 <= `x`, `y` <= 1
  Vec3 renderImagePlanePoint(double x, double y, std::mt19937& thread_mt_19937);
  Vec3 castRay(const Ray& r, int recursion_depth,
               std::mt19937& thread_mt_19937);
  // Returns the direct contribution of `point_lights`
  // to the color going from `point` to direction `out_vector`
  //
  // Assumes that `normal` and `out_vector` are normalized
  Vec3 totalPointLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                            const Material& material);
  // Assumes that `normal` and `out_vector` are normalized
  Vec3 pointLightColor(const Vec3& point, const Vec3& normal,
                       const Vec3& out_vector, const Material& material,
                       const PointLight& point_light);
  // Returns an estimate for the contribution of light from sources
  // other than direct point lights to the color going from
  // `point` to direction `out_vector`
  //
  // Assumes `normal` and `out_vector` are normalized
  Vec3 indirectLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                          const Material& material, int depth,
                          std::mt19937& thread_mt_19937);
  // Samples directions based on `sampling_scheme`, casts ray
  // to that direction and estimates the color of
  // light going to `out_vector`
  Vec3 sampleIndirectLight(const Vec3& point, const Vec3& normal,
                           const Vec3& out_vector, const Material& material,
                           int recursion_depth, std::mt19937& thread_mt_19937);
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
