#ifndef RAYCASTER_RAYCASTER_H
#define RAYCASTER_RAYCASTER_H
#include <random>
#include <vector>

#include "geometry.h"
#include "kd_tree.h"

const double kPi = std::atan(1) * 4;
enum NormalType {
  kSmooth,
  kRough,
};
struct Material {
  Vec3 diffuse;
  Vec3 emitted;
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
class Image {
 public:
  Image(int x_resolution, int y_resolution);
  void setColor(int x_pos, int y_pos, const Vec3& color);
  // scales values to [0, 1]
  void scaleMaxTo(double new_max);
  bool savePPM(const std::string& file);

 private:
  double maxColorValue();
  int bufferPos(int x, int y);
  char floatToByte(double val);
  int x_resolution;
  int y_resolution;
  std::vector<Vec3> buffer;
};
inline void Image::setColor(int x, int y, const Vec3& color) {
  assert(x >= 0 && x < x_resolution);
  assert(y >= 0 && y < y_resolution);
  assert(color[0] > -EPS && color[1] > -EPS && color[2] > -EPS);
  buffer[(size_t)bufferPos(x, y)] = color;
}
inline int Image::bufferPos(int x, int y) { return y * x_resolution + x; }
inline char Image::floatToByte(double val) { return val * 255; }
class Scene {
 public:
  Scene(const Camera& c);
  bool addModelFromFile(const std::string& file, const Vec3& position,
                        NormalType normal_type);
  void addPointLight(const Vec3& position, const Vec3& color);
  Image render(int x_resolution, int y_resolution, int n_rays_per_pixel);

 private:
  Vec3 renderPixel(int x_pixel, int y_pixel, int x_resolution, int y_resolution,
                   int n_rays_per_pixel);
  Vec3 renderImagePlanePoint(double x, double y);
  Vec3 castRay(const Ray& r);
  Model model;
  std::vector<PointLight> point_lights;
  Camera camera;
  Tree kd_tree;
  // used for anti-aliasing
  std::mt19937 mt;
  std::uniform_real_distribution<double> subpixel_sample_distribution;
};
#endif
