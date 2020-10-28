#ifndef RAYCASTER_INPUT_H
#define RAYCASTER_INPUT_H
#include <string>

#include "geometry.h"
#include "raycaster.h"

struct CameraConfig {
  Vec3 pos{0.0};
  Vec3 front{0.0, 0.0, -1.0};
  Vec3 up{0.0, 1.0, 0.0};
  double x_fov = kPi / 2.0;
  double y_fov = kPi / 2.0;
};
struct ModelConfig {
  std::string file;
  Vec3 pos{0.0};
  NormalType normal = NormalType::kRough;
};
struct PointLightConfig {
  Vec3 pos{0.0};
  Vec3 color{1.0};
};
struct EnvironmentLightConfig {
  Vec3 color{0.0};
  EnvironmentLightType type = EnvironmentLightType::kUniform;
  Vec3 direction{0.0, 1.0, 0.0};
  double exp = 1.0;
};
struct RenderConfig {
  int width = 500;
  int height = 500;
  int pixel_rays = 1;
  int depth = 1;
  int branching = 1;
};
struct ImageConfig {
  std::string file;
  double truncate = 1.0;
  double scale_max = 1.0;
};
struct Recipe {
  CameraConfig camera;
  Vector<ModelConfig> models;
  Vector<PointLightConfig> point_lights;
  EnvironmentLightConfig environment_light;
  RenderConfig render;
  Vector<ImageConfig> images;
};
std::ostream& operator<<(std::ostream& out, const CameraConfig& a);
std::ostream& operator<<(std::ostream& out, const ModelConfig& a);
std::ostream& operator<<(std::ostream& out, const PointLightConfig& a);
std::ostream& operator<<(std::ostream& out, const EnvironmentLightConfig& a);
std::ostream& operator<<(std::ostream& out, const RenderConfig& a);
std::ostream& operator<<(std::ostream& out, const ImageConfig& a);
std::ostream& operator<<(std::ostream& out, const Recipe& a);
Recipe loadRecipe(const std::string& file, bool verbose);
#endif
