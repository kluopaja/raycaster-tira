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
  // TODO think if this should be Triangle
  // change after we see where exactly this will be used
  Vec3 normals[3];
  Material material;
};
class Model {
 public:
  Model() = default;
  Model(Model&& a) = default;
  Model& operator=(Model&& a) = default;
  void concatenate(const Model& a);
  void translate(const Vec3& v);
  Tree buildKdTree(double k_t, double k_i);
  std::vector<SceneTriangle> scene_triangles;
  std::vector<Material> materials;
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
// For parallel rendering
struct ThreadResources {
  std::mt19937 thread_mt;
  std::vector<std::pair<int, int> > pixels;
};
class Scene {
 public:
  Scene(const Camera& c);
  bool addModelFromFile(const std::string& file, const Vec3& position,
                        NormalType normal_type);
  void addPointLight(const Vec3& position, const Vec3& color);
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
  Vec3 environment_light_color;
  Model model;
  std::vector<PointLight> point_lights;
  Camera camera;
  Tree kd_tree;
  std::mt19937 mt_19937;
  std::uniform_real_distribution<double> subpixel_sample_distribution;
  int max_recursion_depth;
  int n_recursion_rays;
  SamplingScheme sampling_scheme;
};
#endif
