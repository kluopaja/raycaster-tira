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
#include <tuple>

#include "kd_tree.h"
#include "lighting_utils.h"
#include "model_loader.h"
#include "utils.h"

const int kNThreads = 8;
Vec3 SceneTriangle::normalAt(const Vec2& bary_coords) {
  Vec3 n = Triangle(normals[0], normals[1], normals[2]).pointFromBary(bary_coords);
  n /= n.norm();
  return n;
}
void Model::translate(const Vec3& v) {
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    scene_triangles[i].triangle.translate(v);
  }
}
void Model::concatenate(const Model& a) {
  scene_triangles.insert(scene_triangles.end(), a.scene_triangles.begin(),
                         a.scene_triangles.end());
  materials.insert(materials.end(), a.materials.begin(), a.materials.end());
}
Tree Model::buildKdTree(double k_t, double k_i) {
  std::vector<SceneTriangle*> scene_triangle_pointers;
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    scene_triangle_pointers.push_back(&scene_triangles[i]);
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
Scene::Scene(const Camera& c)
    : environment_light_color(Vec3(0.0)),
      camera(c),
      sampling_scheme(kImportanceSampling) {}
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
void Scene::setEnvironmentLightColor(const Vec3& color) {
  assert(color[0] > -EPS && color[1] > -EPS && color[2] > -EPS);
  environment_light_color = color;
}
void Scene::setSamplingScheme(SamplingScheme s) { sampling_scheme = s; }
Image Scene::render(int x_resolution, int y_resolution, int n_rays_per_pixel,
                    int max_recursion_depth, int n_recursion_rays) {
  kd_tree = model.buildKdTree(1.0, 5.0);
  std::cout << "start rendering... " << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  mt_19937 = std::mt19937(1337);
  subpixel_sample_distribution = std::uniform_real_distribution(0.0, 1.0);
  this->max_recursion_depth = max_recursion_depth;
  this->n_recursion_rays = n_recursion_rays;
#if ENABLE_PARALLEL
  // constructing a mt19937 seems to be quite expensive so only construct
  // kNThreads of them
  std::vector<ThreadResources> thread_resources;
  for (int i = 0; i < kNThreads; ++i) {
    unsigned int seed = mt_19937();
    thread_resources.push_back({std::mt19937(seed), {}});
  }
  for (int x_pixel = 0; x_pixel < x_resolution; ++x_pixel) {
    for (int y_pixel = 0; y_pixel < y_resolution; ++y_pixel) {
      int target_thread = (x_pixel + y_pixel * y_resolution) % kNThreads;
      thread_resources[(size_t)target_thread].pixels.push_back(
          {x_pixel, y_pixel});
    }
  }
  Image image(x_resolution, y_resolution);
  for_each(std::execution::par, thread_resources.begin(),
           thread_resources.end(), [&](ThreadResources& x) {
             for (auto pixel : x.pixels) {
               Vec3 pixel_color =
                   renderPixel(pixel.first, pixel.second, x_resolution,
                               y_resolution, n_rays_per_pixel, x.thread_mt);
               image.setColor(pixel.first, pixel.second, pixel_color);
             }
           });
#else
  Image image(x_resolution, y_resolution);
  for (int x_pixel = 0; x_pixel < x_resolution; ++x_pixel) {
    for (int y_pixel = 0; y_pixel < y_resolution; ++y_pixel) {
      // TODO make this equivalent with the parallel part
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
    return environmentLightColor(r);
  }
  Vec3 intersection_point =
      sp.scene_triangle->triangle.pointFromBary(sp.bary_coords);
  Vec3 normal = sp.scene_triangle->normalAt(sp.bary_coords);

  Vec3 light_color(0.0);
  light_color = light_color + sp.scene_triangle->material.emitted;
  // point lights
  light_color = light_color + totalPointLightColor(intersection_point, normal,
                                                   -1.0 * r.direction,
                                                   sp.scene_triangle->material);
  if (recursion_depth < max_recursion_depth) {
    light_color =
        light_color + indirectLightColor(intersection_point, normal,
                                         -1.0 * r.direction,
                                         sp.scene_triangle->material,
                                         recursion_depth, thread_mt_19937);
  }
  assert(light_color[0] > -EPS && light_color[1] > -EPS &&
         light_color[2] > -EPS);
  return light_color;
}
Vec3 Scene::totalPointLightColor(const Vec3& point, Vec3 normal,
                                 Vec3 out_vector, const Material& material) {
  normal /= normal.norm();
  out_vector /= out_vector.norm();
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
  Vec3 shadow_ray_direction = point_light.position - point;
  Vec3 shadow_ray_start = point + shadow_ray_direction * 100 * EPS;
  if (kd_tree.trianglesIntersectSegment(shadow_ray_start,
                                        point_light.position)) {
    return Vec3(0.0);
  }
  Vec3 light_direction = point_light.position - point;
  double distance_dimming = 1 / light_direction.dot(light_direction);
  Vec3 in_vector = light_direction / light_direction.norm();
  Vec3 in_color = distance_dimming * point_light.color;
  Vec3 out_color = material.apply(in_vector, normal, out_vector, in_color);
  if (out_color.norm() > 10) {
    std::cout << "distance dimmed point: " << in_color << std::endl;
    std::cout << out_color << std::endl;
  }
  return material.apply(in_vector, normal, out_vector, in_color);
}
Vec3 Scene::indirectLightColor(const Vec3& point, Vec3 normal, Vec3 out_vector,
                               const Material& material, int recursion_depth,
                               std::mt19937& thread_mt_19937) {
  normal /= normal.norm();
  out_vector /= out_vector.norm();
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
    in_vector = uniformRandomSpherePoint(normal, thread_mt_19937);
    // std::cerr << in_vector << normal << out_vector << std::endl;
    in_pdf = 1 / (4 * kPi);
  }
  assert(in_pdf > EPS);

  Ray in_ray(point + in_vector * 100 * EPS, in_vector);
  Vec3 in_color = castRay(in_ray, recursion_depth + 1, thread_mt_19937);
  return material.apply(in_vector, normal, out_vector, in_color) / in_pdf;
}
Vec3 Scene::environmentLightColor(const Ray& r) {
  // "unused parameter r warning" can be ignored for now
  return environment_light_color;
}
