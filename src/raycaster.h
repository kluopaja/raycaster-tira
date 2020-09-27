#ifndef RAYCASTER_RAYCASTER_H
#define RAYCASTER_RAYCASTER_H
#include "geometry.h"
struct Material {
  Vec3 diffuse;
  Vec3 specular;
  Vec3 emitted;
};
struct SceneTriangle {
  Triangle triangle;
  Vec3 normals[3];
  Material* material;
};
#endif
