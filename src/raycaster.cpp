#include "raycaster.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>

#include "kd_tree.h"
#include "model_loader.h"

Model::Model(Model&& a) noexcept
    : scene_triangles(std::move(a.scene_triangles)),
      materials(std::move(a.materials)) {
  std::cout << "move constructing Model..." << std::endl;
}
Model& Model::operator=(Model&& a) noexcept {
  // detect self-assignment
  if (&a == this) {
    return *this;
  }
  for (auto x : scene_triangles) {
    delete x;
  }
  for (auto x : materials) {
    delete x;
  }
  scene_triangles = std::move(a.scene_triangles);
  materials = std::move(a.materials);
  return *this;
}
void Model::translate(const Vec3& v) {
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    scene_triangles[i]->triangle.translate(v);
  }
}
Model::~Model() {
  for (auto x : scene_triangles) {
    delete x;
  }
  for (auto x : materials) {
    delete x;
  }
}
void Model::concatenate(Model&& a) {
  // detect self-concatenation
  if (&a == this) {
    return;
  }
  scene_triangles.insert(scene_triangles.end(), a.scene_triangles.begin(),
                         a.scene_triangles.end());
  materials.insert(materials.end(), a.materials.begin(), a.materials.end());

  // this will prevent deleting SceneTriangles or Materials
  // when destructor of 'a' is called
  a.scene_triangles.clear();
  a.materials.clear();
}
PointLight::PointLight(const Vec3& position, const Vec3& color)
    : position(position), color(color) {
  assert(color[0] > 0 && color[1] > 0 && color[2] > 0);
}
Camera::Camera(const Vec3& position, const Vec3& front_vec, const Vec3& up_vec,
               double x_fov, double y_fov)
    : position(position) {
  assert(EPS < x_fov && x_fov < kPi - EPS);
  assert(EPS < y_fov && y_fov < kPi - EPS);
  assert(front_vec.norm() > EPS);
  assert(up_vec.norm() > EPS);
  assert(std::abs(front_vec.dot(up_vec)) < EPS);

  image_down_vec =
      -2 * up_vec / up_vec.norm() * std::tan(y_fov / 2) * front_vec.norm();
  Vec3 right_vec = front_vec.cross(up_vec);
  image_right_vec =
      2 * right_vec / right_vec.norm() * std::tan(x_fov / 2) * front_vec.norm();
  image_position =
      position + front_vec - image_down_vec / 2 - image_right_vec / 2;
}
Ray Camera::rayFromImagePlane(double x, double y) {
  assert(x > -EPS && x < 1 + EPS);
  assert(y > -EPS && y < 1 + EPS);
  Vec3 ray_start = image_position + x * image_right_vec + y * image_down_vec;
  return Ray(ray_start, ray_start - position);
}
Scene::Scene(const Camera& c) : camera(c) {}
bool Scene::addModelFromFile(const std::string& file, const Vec3& position,
                             NormalType normal_type) {
  Model new_model;
  if (!loadModel(file, new_model, normal_type)) {
    return false;
  }
  new_model.translate(position);
  model.concatenate(std::move(new_model));
  return true;
}
void Scene::addPointLight(const Vec3& position, const Vec3& color) {
  point_lights.emplace_back(position, color);
}
Image Scene::render(int x_resolution, int y_resolution, int n_rays_per_pixel) {
  kd_tree = buildKdTree(model.scene_triangles, 1.0, 5.0);
  std::cerr << "start rendering... " << std::endl;
  mt = std::mt19937(1337);
  subpixel_sample_distribution = std::uniform_real_distribution(0.0, 1.0);
  Image image(x_resolution, y_resolution);
  for (int x_pixel = 0; x_pixel < x_resolution; ++x_pixel) {
    for (int y_pixel = 0; y_pixel < y_resolution; ++y_pixel) {
      Vec3 pixel_color = renderPixel(x_pixel, y_pixel, x_resolution,
                                     y_resolution, n_rays_per_pixel);
      for (int i = 0; i < n_rays_per_pixel; ++i) {
        // double x = (x_pixel + subpixel_ray_distribution(mt)) / x_resolution;
        // double y = (y_pixel + subpixel_ray_distribution(mt)) / y_resolution;
        double x = (x_pixel + 0.5) / x_resolution;
        double y = (y_pixel + 0.5) / y_resolution;
        Ray r = camera.rayFromImagePlane(x, y);
        pixel_color = pixel_color + castRay(r);
      }
      pixel_color = pixel_color / n_rays_per_pixel;
      image.setColor(x_pixel, y_pixel, pixel_color);
    }
  }
  return image;
}
Vec3 Scene::renderPixel(int x_pixel, int y_pixel, int x_resolution,
                        int y_resolution, int n_rays_per_pixel) {
  assert(n_rays_per_pixel > 0);
  assert(x_resolution > 0);
  assert(y_resolution > 0);
  if (n_rays_per_pixel == 1) {
    // use deterministic sampling
    double x = (x_pixel + 0.5) / x_resolution;
    double y = (y_pixel + 0.5) / y_resolution;
    return renderImagePlanePoint(x, y);
  } else {
    Vec3 color(0.0);
    for (int i = 0; i < n_rays_per_pixel; ++i) {
      double x = (x_pixel + subpixel_sample_distribution(mt)) / x_resolution;
      double y = (y_pixel + subpixel_sample_distribution(mt)) / y_resolution;
      color = color + renderImagePlanePoint(x, y);
    }
    return color;
  }
}
Vec3 Scene::renderImagePlanePoint(double x, double y) {
  Ray r = camera.rayFromImagePlane(x, y);
  return castRay(r);
}
Vec3 Scene::castRay(const Ray& r) {
  ScenePoint sp = kd_tree.getClosestRayIntersection(r);
  // no intersection
  if (sp.scene_triangle == nullptr) {
    return Vec3(0.0, 0.0, 0.0);
  }
  Vec3 intersection_point =
      sp.scene_triangle->triangle.pointFromBary(sp.bary_coords);
  // TODO change to something else
  Vec3 normal =
      Triangle(sp.scene_triangle->normals[0], sp.scene_triangle->normals[1],
               sp.scene_triangle->normals[2])
          .pointFromBary(sp.bary_coords);
  normal = normal / normal.norm();

  // std::cerr << normal.norm() << std::endl;
  assert(std::abs(normal.norm() - 1.0) < EPS);
  Vec3 light_color(0.0);
  // std::cerr << "found triangle! " << std::endl;
  // std::cerr << "ray: " << r << std::endl;
  // std::cerr << "intersection_point " << intersection_point << std::endl;
  // point lights
  for (size_t i = 0; i < point_lights.size(); ++i) {
    Vec3 shadow_ray_direction = point_lights[i].position - intersection_point;
    Vec3 shadow_ray_start =
        intersection_point + shadow_ray_direction * 100 * EPS;
    // std::cerr << "shadow ray: " << shadow_ray_start
    //           << " " << shadow_ray_direction << std::endl;
    if (kd_tree.trianglesIntersectSegment(shadow_ray_start,
                                          point_lights[i].position)) {
      continue;
    }
    // std::cerr << "found light! " << std::endl;
    Vec3 light_vector = point_lights[i].position - intersection_point;
    // std::cerr << "light_vector: " << light_vector << std::endl;
    // std::cerr << "normal: " << normal << std::endl;
    double light_distance = light_vector.norm();
    light_vector = light_vector / light_distance;
    // diffuse light
    // no need to use any pi term here to take care of the conservation of
    // energy because the light intensity
    // is anyway just some random constant
    light_color = light_color + std::max(0.0, normal.dot(light_vector)) /
                                    light_distance *
                                    sp.scene_triangle->material->diffuse;
  }
  return light_color;
}
Image::Image(int x_resolution, int y_resolution)
    : x_resolution(x_resolution), y_resolution(y_resolution) {
  assert(x_resolution > 0);
  assert(y_resolution > 0);
  buffer.resize((size_t)x_resolution * (size_t)y_resolution);
}
void Image::scaleMaxTo(double new_max) {
  assert(new_max > EPS);
  double old_max = maxColorValue();
  if (old_max < EPS) {
    return;
  }
  double scale_factor = new_max / old_max;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      size_t pos = (size_t)bufferPos(x, y);
      buffer[pos] = buffer[pos] * scale_factor;
    }
  }
}
bool Image::savePPM(const std::string& file) {
  std::ofstream fout(file, std::ios::binary);
  if (!fout.is_open()) {
    return false;
  }
  // binary mode
  fout << "P6" << std::endl;
  fout << x_resolution << " " << y_resolution << std::endl;
  fout << "255" << std::endl;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      // r g b
      for (int i = 0; i < 3; ++i) {
        fout << floatToByte(buffer[(size_t)bufferPos(x, y)][i]);
      }
    }
  }
  fout.close();
  return true;
}
double Image::maxColorValue() {
  double max_value = 0;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      for (int i = 0; i < 3; ++i) {
        max_value = std::max(max_value, buffer[(size_t)bufferPos(x, y)][i]);
      }
    }
  }
  return max_value;
}
