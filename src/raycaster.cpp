#ifndef ENABLE_PARALLEL
#define ENABLE_PARALLEL 0
#endif
#include "raycaster.h"

#include <cassert>
#include <chrono>
#include <cmath>

#if ENABLE_PARALLEL
#include <execution>
#endif

#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <tuple>

#include "kd_tree.h"
#include "lighting_utils.h"
#include "model_loader.h"
#include "utils.h"

const int kNThreads = 8;
Vec3 SceneTriangle::normalAt(const Vec2& bary_coords) {
  Vec3 n =
      Triangle(normals[0], normals[1], normals[2]).pointFromBary(bary_coords);
  n /= n.norm();
  return n;
}
void Model::translate(const Vec3& v) {
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    scene_triangles[i].triangle.translate(v);
  }
}
void Model::concatenate(const Model& a) {
  for (auto x : a.scene_triangles) {
    scene_triangles.pushBack(x);
  }
  for (auto x : a.materials) {
    materials.pushBack(x);
  }
}
Tree Model::buildKdTree(double k_t, double k_i) {
  Vector<SceneTriangle*> scene_triangle_pointers;
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    scene_triangle_pointers.pushBack(&scene_triangles[i]);
  }
  return ::buildKdTree(scene_triangle_pointers, k_t, k_i);
}
PointLight::PointLight(const Vec3& position, const Vec3& color)
    : position(position), color(color) {
  assert(color[0] > -EPS && color[1] > -EPS && color[2] > -EPS);
}
Camera::Camera(const Vec3& position, const Vec3& front_vec, const Vec3& up_vec,
               double x_fov, double y_fov)
    : position(position) {
  assert(EPS < x_fov && x_fov < kPi - EPS);
  assert(EPS < y_fov && y_fov < kPi - EPS);
  assert(front_vec.norm() > EPS);
  assert(up_vec.norm() > EPS);
  assert(std::abs(front_vec.dot(up_vec)) < EPS &&
         "Camera up and front vectors should be perpendicular");
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
void EnvironmentLight::setUniform() { is_directed = 0; }
void EnvironmentLight::setDirected(const Vec3& direction, double exponent) {
  assert(direction.norm() > EPS);
  assert(exponent >= 0);
  is_directed = 1;
  this->direction = direction;
  this->direction /= this->direction.norm();
  cosine_exp = exponent;
}
void EnvironmentLight::setColor(const Vec3& color) { this->color = color; }
Scene::Scene(const Camera& c)
    : camera(c), sampling_scheme(kImportanceSampling) {}
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
  point_lights.pushBack(PointLight(position, color));
}
void Scene::setEnvironmentLightUniform() { environment_light.setUniform(); }
void Scene::setEnvironmentLightDirected(const Vec3& direction,
                                        double exponent) {
  assert(direction.norm() > EPS &&
         "Environment light direction should not "
         "be zero");
  assert(exponent >= 0 && "Exponent should be non-negative");
  environment_light.setDirected(direction, exponent);
}
void Scene::setEnvironmentLightColor(const Vec3& color) {
  assert(color[0] > -EPS && color[1] > -EPS && color[2] > -EPS &&
         "Environment light color should not be negative");
  environment_light.setColor(color);
}
void Scene::setSamplingScheme(SamplingScheme s) { sampling_scheme = s; }
void Scene::setCosts(double new_k_t, double new_k_i) {
  assert(new_k_t > 0);
  assert(new_k_i > 0);
  k_t = new_k_t;
  k_i = new_k_i;
}
Image Scene::render(int x_resolution, int y_resolution, int n_rays_per_pixel,
                    int max_recursion_depth, int n_recursion_rays) {
  assert(x_resolution > 0);
  assert(y_resolution > 0);
  assert(n_rays_per_pixel > 0);
  assert(max_recursion_depth > 0);
  assert(n_recursion_rays > 0);

  std::cout << "starting to build the kd tree..." << std::endl;
  kd_tree = model.buildKdTree(k_t, k_i);

  std::cout << "starting to render... " << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  mt_19937 = std::mt19937(1337);
  subpixel_sample_distribution = std::uniform_real_distribution(0.0, 1.0);

  this->max_recursion_depth = max_recursion_depth;
  this->n_recursion_rays = n_recursion_rays;
#if ENABLE_PARALLEL
  // constructing a mt19937 seems to be quite expensive so only construct
  // kNThreads of them
  Vector<ThreadResources> thread_resources;
  for (int i = 0; i < kNThreads; ++i) {
    unsigned int seed = mt_19937();
    thread_resources.pushBack({std::mt19937(seed), {}});
  }
  for (int x_pixel = 0; x_pixel < x_resolution; ++x_pixel) {
    for (int y_pixel = 0; y_pixel < y_resolution; ++y_pixel) {
      int target_thread = (x_pixel + y_pixel * y_resolution) % kNThreads;
      thread_resources[(size_t)target_thread].pixels.pushBack(
          {x_pixel, y_pixel});
    }
  }
  Image image(x_resolution, y_resolution);
  std::for_each(std::execution::par, thread_resources.begin(),
                thread_resources.end(), [&](ThreadResources& x) {
                  for (auto pixel : x.pixels) {
                    Vec3 pixel_color = renderPixel(
                        pixel.first, pixel.second, x_resolution, y_resolution,
                        n_rays_per_pixel, x.thread_mt);
                    image.setColor(pixel.first, pixel.second, pixel_color);
                  }
                });
#else
  Image image(x_resolution, y_resolution);
  for (int x_pixel = 0; x_pixel < x_resolution; ++x_pixel) {
    for (int y_pixel = 0; y_pixel < y_resolution; ++y_pixel) {
      // Note that this is not exactly equivalent to the parallelel
      // version. (Random numbers can be different)
      Vec3 pixel_color = renderPixel(x_pixel, y_pixel, x_resolution,
                                     y_resolution, n_rays_per_pixel, mt_19937);
      image.setColor(x_pixel, y_pixel, pixel_color);
    }
  }
#endif
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> render_time = end_time - start_time;
  std::cout << "rendering took: " << render_time.count() << " s" << std::endl;
  return image;
}
Vec3 Scene::renderPixel(int x_pixel, int y_pixel, int x_resolution,
                        int y_resolution, int n_rays_per_pixel,
                        std::mt19937& thread_mt_19937) {
  assert(n_rays_per_pixel > 0);
  assert(x_resolution > 0);
  assert(y_resolution > 0);
  if (n_rays_per_pixel == 1) {
    // use deterministic sampling
    double x = (x_pixel + 0.5) / x_resolution;
    double y = (y_pixel + 0.5) / y_resolution;
    return renderImagePlanePoint(x, y, mt_19937);
  } else {
    Vec3 color(0.0);
    for (int i = 0; i < n_rays_per_pixel; ++i) {
      double x = (x_pixel + subpixel_sample_distribution(thread_mt_19937)) /
                 x_resolution;
      double y = (y_pixel + subpixel_sample_distribution(thread_mt_19937)) /
                 y_resolution;
      color = color + renderImagePlanePoint(x, y, thread_mt_19937);
    }
    return color / n_rays_per_pixel;
  }
}
Vec3 Scene::renderImagePlanePoint(double x, double y,
                                  std::mt19937& thread_mt_19937) {
  Ray r = camera.rayFromImagePlane(x, y);
  return castRay(r, 1, thread_mt_19937);
}
Vec3 Scene::castRay(const Ray& r, int recursion_depth,
                    std::mt19937& thread_mt_19937) {
  ScenePoint sp = kd_tree.getClosestRayIntersection(r);
  // no intersection
  if (sp.scene_triangle == nullptr) {
    return environment_light.colorAtDirection(r.getDirection());
  }
  Vec3 intersection_point =
      sp.scene_triangle->triangle.pointFromBary(sp.bary_coords);
  Vec3 normal = sp.scene_triangle->normalAt(sp.bary_coords);
  Vec3 out_vector = normalize(-1.0 * r.getDirection());

  Vec3 light_color(0.0);
  light_color = light_color + sp.scene_triangle->material.emitted;
  // point lights
  light_color =
      light_color + totalPointLightColor(intersection_point, normal, out_vector,
                                         sp.scene_triangle->material);
  if (recursion_depth < max_recursion_depth) {
    light_color =
        light_color + indirectLightColor(intersection_point, normal, out_vector,
                                         sp.scene_triangle->material,
                                         recursion_depth, thread_mt_19937);
  }
  assert(light_color[0] > -EPS && light_color[1] > -EPS &&
         light_color[2] > -EPS);
  return light_color;
}
Vec3 Scene::totalPointLightColor(const Vec3& point, Vec3 normal,
                                 Vec3 out_vector, const Material& material) {
  Vec3 point_color(0.0);
  for (size_t i = 0; i < point_lights.size(); ++i) {
    point_color = point_color + pointLightColor(point, normal, out_vector,
                                                material, point_lights[i]);
  }
  return point_color;
}
// normal and out_vector should be normalized
Vec3 Scene::pointLightColor(const Vec3& point, const Vec3& normal,
                            const Vec3& out_vector, const Material& material,
                            const PointLight& point_light) {
  assert(std::abs(normal.norm() - 1.0) < EPS);
  assert(std::abs(out_vector.norm() - 1.0) < EPS);
  Vec3 light_direction = point_light.position - point;
  Vec3 in_vector = normalize(light_direction);
  // To prevent self-collisions, add some offset before
  // the start of the ray
  Vec3 shadow_ray_start = point + in_vector * 100.0 * EPS;
  if (kd_tree.trianglesIntersectSegment(shadow_ray_start,
                                        point_light.position)) {
    return Vec3(0.0);
  }
  double distance_dimming = 1.0 / light_direction.dot(light_direction);
  Vec3 in_color = distance_dimming * point_light.color;
  return material.apply(in_vector, normal, out_vector, in_color);
}
Vec3 Scene::indirectLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                               const Material& material, int recursion_depth,
                               std::mt19937& thread_mt_19937) {
  assert(std::abs(normal.norm() - 1.0) < EPS);
  assert(std::abs(out_vector.norm() - 1.0) < EPS);
  Vec3 out_color(0.0);
  for (int i = 0; i < n_recursion_rays; ++i) {
    out_color =
        out_color + sampleIndirectLight(point, normal, out_vector, material,
                                        recursion_depth, thread_mt_19937);
  }
  // calculate final integral estimate
  return out_color / n_recursion_rays;
}
Vec3 Scene::sampleIndirectLight(const Vec3& point, const Vec3& normal,
                                const Vec3& out_vector,
                                const Material& material, int recursion_depth,
                                std::mt19937& thread_mt_19937) {
  Vec3 in_vector;
  double in_pdf;
  if (sampling_scheme == kImportanceSampling) {
    in_vector = material.importanceSample(normal, out_vector, thread_mt_19937);
    in_pdf = material.importanceSamplePdf(in_vector, normal, out_vector);
  } else {
    in_vector = uniformRandomSpherePoint(thread_mt_19937);
    in_pdf = 1 / (4 * kPi);
  }
  assert(in_pdf > EPS);
  Vec3 ray_start = point + 100 * EPS * in_vector;
  ;

  Ray in_ray(ray_start, in_vector);
  Vec3 in_color = castRay(in_ray, recursion_depth + 1, thread_mt_19937);
  return material.apply(in_vector, normal, out_vector, in_color) / in_pdf;
}
