#include "raycaster.h"
#include <memory>
#include <iostream>

/*
class Camera {
  Camera(const Vec3& position, const Vec3& direction, const Vec3& fov,
         const Vec2& resolution): position(position), direction(direction),
                                  fov(fov), resolution(resolution);
  Vec3 position;
  Vec3 direction;
  Vec2 fov; // fied of view in degrees
  Vec2 resolution;
}
*/
Model::Model(Model&& a) noexcept: scene_triangles(std::move(a.scene_triangles)),
                         materials(std::move(a.materials)) {
  std::cout << "move constructing Model..." << std::endl;
}
Model& Model::operator=(Model&& a) noexcept {
  for(auto x: scene_triangles) {
    delete x;
  }
  for(auto x: materials) {
    delete x;
  }
  scene_triangles = std::move(a.scene_triangles);
  materials = std::move(a.materials);
  return *this;
}
Model::~Model() {
  for(auto x: scene_triangles) {
    delete x;
  }
  for(auto x: materials) {
    delete x;
  }
}

