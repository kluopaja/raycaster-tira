#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <string>
#include <random>

#include "../../utils.h"
#include "../../vector.h"
#include "performance_utils.h"
#include "../test_utils.h"
#include "../../raycaster.h"
#include "../../model_loader.h"

namespace {
void generateKdStructureReport() {
  Model bunny_model;
  if(!loadModel("../models/bunny/bunny_hires.obj", bunny_model,
                NormalType::kRough)) {
    std::cerr << "Failed loading bunny model" << std::endl;
    std::exit(1);
  }
  std::cout << "Generating a kd tree structure report." << std::endl;
  
  std::cout << "From bunny_hires.obj with"
               " parameters k_t = 15, k_i = 20..." << std::endl;
  Timer timer;
  timer.start();
  Tree bunny_tree = bunny_model.buildKdTree(15.0, 20.0);
  double time = timer.elapsed();
  std::cout << "time: " << time/1e6 << std::endl;
  std::cout << "numLeaves: " << bunny_tree.numLeaves() << std::endl;
  std::cout << "triangles: " << bunny_model.scene_triangles.size() << std::endl;
  std::cout << "nonEmptyLeaves: " << bunny_tree.numNonEmptyLeaves() << std::endl;
  std::cout << "cost " << bunny_tree.cost(15, 20) << std::endl;
  std::cout << "average triangles in non-empty leaves: " << bunny_tree.averageTrianglesInLeaf() << std::endl;
}
struct RunTimes {
  double kd_tree_time = 0;
  double naive_time = 0;
};
RunTimes operator+(const RunTimes& a, const RunTimes& b) {
  return {a.kd_tree_time + b.kd_tree_time, a.naive_time + b.naive_time};
}
RunTimes operator/(const RunTimes& a, double b) {
  return {a.kd_tree_time / b, a.naive_time / b};
}
std::ostream& operator<<(std::ostream& out, const RunTimes& a) {
  out << std::left;
  out << std::setw(25);
  out << "Kd tree time: " << a.kd_tree_time/1e6 << "ms" << std::endl;
  out << std::setw(25);
  out << "naive time: " << a.naive_time/1e6 << "ms" << std::endl;
  return out;
}
RunTimes runRandomQueries(const std::string& file, int n_queries) {
  Model model;
  if(!loadModel(file, model, NormalType::kRough)) {
    std::cerr << "Failed loading model" << std::endl;
    std::exit(1);
  }
  Tree tree = model.buildKdTree(15, 20);
  Voxel bb = tree.boundingBox();
  std::mt19937 mt(1337);
  Vector<Ray> rays;
  for(int i = 0; i < n_queries; ++i) {
    Vec3 ray_start;
    for(int j = 0; j < 3; ++j) {
      std::uniform_real_distribution dist(bb.lo[j], bb.hi[j]);
      ray_start[j] = dist(mt);
    }
    Ray r(ray_start,
          uniformRandomSpherePoint(mt)/100);
    rays.pushBack(r);
  }

  // Extract Triangles for firstRayTriangleIntersection
  Vector<SceneTriangle*> scene_p;
  for(auto& x : model.scene_triangles) {
    scene_p.pushBack(&x);
  }
  Vector<Triangle> triangles = extractTriangles(scene_p);
  RunTimes query_times;
  std::cout << "triangles: " << triangles.size() << std::endl;
  Timer timer;
  timer.start();
  int n_intersections_naive = 0;
  for(auto ray: rays) {
    RayTriangleIntersection rti = firstRayTriangleIntersection(triangles, ray);
    if(rti.index < triangles.size()) {
      ++n_intersections_naive;
    }
  }
  query_times.naive_time = timer.elapsed();
  timer.start();
  int n_intersections_kd_tree = 0;
  for(auto ray: rays) {
    ScenePoint sp = tree.getClosestRayIntersection(ray);
    if(sp.scene_triangle != nullptr) {
      ++n_intersections_kd_tree;
    }
  }
  query_times.kd_tree_time = timer.elapsed();
  assert(n_intersections_naive == n_intersections_kd_tree);
  std::cout << query_times << std::endl;
  return query_times;
}
void generateKdQueryReport() {

  std::cout << "**********************" << std::endl;
  std::cout << "Generating Kd tree query report..." << std::endl;
  RunTimes query_times = {0.0, 0.0};
  for(int i = 0; i < 10; ++i) {
    RunTimes tmp = runRandomQueries("../models/armadillo/armadillo.obj",
                                            100);
    query_times = query_times + tmp;
  }
  query_times = query_times / 10;
  std::cout << "100 queries for armadillo.obj (mean of 10 repeats):\n";
  std::cout << query_times << std::endl;
}
void generateKdRenderReport() {
  std::cout << "**********************" << std::endl;
  std::cout << "Generating Kd tree rendering report..." << std::endl;
  Camera camera(Vec3(0.0, 0.1, 6.0), Vec3(0.0, 0.0, -0.1),
                Vec3(0.0, 1.0, 0.0), kPi/4.0, kPi/4.0);
  Scene scene(camera);
  scene.addModelFromFile("../models/bunny/bunny_hires.obj", 
                         Vec3(0.0, 0.0, 0.0), NormalType::kRough);
  scene.addPointLight(Vec3(5.0, 5.0, 5.0), Vec3(20.0, 20.0, 20.0));
  Timer timer;
  RunTimes render_times;
  timer.start();
  Image result1 = scene.render(100, 100, 1, 6, 2);
  render_times.kd_tree_time = timer.elapsed();
  // Set tree traversal cost to very high to generate
  // only one node to the tree
  scene.setCosts(1e20, 2*EPS);
  timer.start();
  Image result2 = scene.render(100, 100, 1, 6, 2);
  render_times.naive_time = timer.elapsed();
  // check that the results are the same
  assert(result1.distanceTo(result2) < EPS);
  std::cout << "total rendering times (including building time of the kd tree): " << std::endl;
  std::cout << render_times << std::endl;
}
}  // namespace
// Prints performance test results to std::cout
void generateKdReport() {
  std::cout << std::string("*", 20) << std::endl;
  std::cout << "Kd tree performance report" << std::endl;

  generateKdStructureReport();
  generateKdQueryReport();
  generateKdRenderReport();
}
