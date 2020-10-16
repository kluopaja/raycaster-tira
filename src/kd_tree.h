#ifndef RAYCASTER_KD_TREE_H
#define RAYCASTER_KD_TREE_H
#include <memory>
#include <vector>

#include "geometry.h"
// forward declare SceneTriangle to solve circular include issue
struct SceneTriangle;
struct ScenePoint {
  SceneTriangle* scene_triangle;
  Vec2 bary_coords;
};
class Node {
 public:
  Node(std::unique_ptr<Node> left, std::unique_ptr<Node> right,
       const Voxel& voxel, const AxisPlane& plane);
  Node(const Vector<Triangle>& triangles,
       const Vector<SceneTriangle*> scene_triangles, const Voxel& voxel);
  ScenePoint getClosestRayIntersection(const Ray& r) const;

 private:
  // the side with smaller coordinates
  // for example if the plane axis = 0,
  // then maximum x coordinate of left.voxel is at most the
  // minimum x coordinate of right.voxel
  std::unique_ptr<Node> left;
  std::unique_ptr<Node> right;
  Vector<Triangle> triangles;
  Vector<SceneTriangle*> scene_triangles;
  Voxel voxel;
  AxisPlane plane;
  bool isLeaf() const;
};
inline bool Node::isLeaf() const { return (!left && !right); }
class Tree {
 public:
  ScenePoint getClosestRayIntersection(const Ray& r) const;
  bool trianglesIntersectSegment(const Vec3& a, const Vec3& b) const;
  Tree(std::unique_ptr<Node> root) : root(std::move(root)) {}
  Tree() = default;
  // move constructor
  Tree(Tree&& a) noexcept;
  Tree& operator=(Tree&& a) noexcept;
  Tree(const Tree&) = delete;
  Tree& operator=(const Tree&) = delete;

 private:
  std::unique_ptr<Node> root;
};

Tree buildKdTree(const Vector<SceneTriangle*>& scene_triangles, double k_t,
                 double k_i);
Vector<Triangle> extractTriangles(
    const Vector<SceneTriangle*>& scene_triangles);

// Calculates the surface area heuristic cost
// for splits [(n_left + n_plane), (n_right)]
// and [(n_left), (n_right + n_plane)]

// returns the cost of the better split and
// the side to which n_plane was assigned
// in this split
std::pair<double, bool> surfaceAreaHeuristic(double l_area, double r_area,
                                             int n_left, int n_plane,
                                             int n_right, double traversal_cost,
                                             double intersection_cost);
std::pair<double, double> relativeSubvoxelAreas(const Voxel& voxel,
                                                const AxisPlane& plane);
#endif
