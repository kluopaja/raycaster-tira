#ifndef RAYCASTER_RAYCASTER_H
#define RAYCASTER_RAYCASTER_H
#include "geometry.h"
struct Material {
  Vec3 diffuse;
  Vec3 emitted;
};
struct SceneTriangle {
  Triangle triangle;
  Vec3 normals[3];
  Material* material;
};
class Model {
 public:
  Model() = default;
  Model(Model&& a) noexcept;
  Model& operator=(Model&& a) noexcept;
  ~Model();
  std::vector<SceneTriangle*> scene_triangles;
  // necessary for freeing the memory
  std::vector<Material*> materials;
};
// struct Scene {
//   Model model;
// };
#endif
