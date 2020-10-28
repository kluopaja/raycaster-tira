#ifndef RAYCASTER_KD_TREE_H
#define RAYCASTER_KD_TREE_H
#include <memory>

#include "geometry.h"
#include "vector.h"
// forward declare SceneTriangle to solve circular include issue
struct SceneTriangle;
// A struct returned for returning the intersection point from kd tree
struct ScenePoint {
  SceneTriangle* scene_triangle;
  // Barycentric coordinates for the intersection
  Vec2 bary_coords;
};
// A class representing a kd tree node
// Should always be accessed through the Tree class
class Node {
 public:
  Node(std::unique_ptr<Node> left, std::unique_ptr<Node> right,
       const Voxel& voxel, const AxisPlane& plane);
  Node(const Vector<Triangle>& triangles,
       const Vector<SceneTriangle*>& scene_triangles, const Voxel& voxel);
  // Returns the information about the closest triangle in the
  // subtree *this
  //
  // Return {nullptr, Vec2(0.0, 0.0)} if no intersection was found
  ScenePoint getClosestRayIntersection(const Ray& r) const;
  // Utilities used to analyze the tree structure
  //
  // Number of leaf nodes in the subtree
  std::size_t numLeaves() const;
  // Number of non-empty leaf nodes in the subtree
  std::size_t numNonEmptyLeaves() const;
  // Total number of triangles in the subtree
  std::size_t totalTriangles() const;
  // Cost of the subtree according to the equation 3 of the paper
  double cost(double traversal_cost, double intersection_cost) const;
  // Returns the smallest axis aligned bounding box of *this
  Voxel boundingBox() const;

 private:
  // the side with smaller coordinates
  // for example if the plane axis = 0,
  // then maximum x coordinate of left.voxel is at most the
  // minimum x coordinate of right.voxel
  //
  // Every node of the tree should have exactly 0 or 2 children.
  std::unique_ptr<Node> left;
  std::unique_ptr<Node> right;
  // The triangles are copied to every subtree. This might
  // increase the memory consumption
  // but was done to speed up the code
  //
  // Didn't really seem to have a big effect
  Vector<Triangle> triangles;
  // Pointers to SceneTriangle objects which store
  // information such as materials and normals
  Vector<SceneTriangle*> scene_triangles;
  // The area covered by the subtree
  Voxel voxel;
  // The plane dividing `voxel` into subtree voxels
  AxisPlane plane;
  bool isLeaf() const;
};
inline bool Node::isLeaf() const { return (!left && !right); }
// A class representing a kd tree
class Tree {
 public:
  // Returns information about the first triangle intersected
  // by ray `r`
  ScenePoint getClosestRayIntersection(const Ray& r) const;
  // Returns 1 if segment from `a` to `b` intersects a triangle stored in
  // the kd tree
  bool trianglesIntersectSegment(const Vec3& a, const Vec3& b) const;
  Tree(std::unique_ptr<Node> root) : root(std::move(root)) {}
  Tree() = default;
  // Utilities used to analyze the tree structure
  // Number of leaf nodes in the tree
  std::size_t numLeaves() const;
  // Number of non-empty leaf nodes in the tree
  std::size_t numNonEmptyLeaves() const;
  // Average number of triangles in non-empty leaves
  double averageTrianglesInLeaf() const;
  // Cost of the tree according to the equation 3 of the paper
  double cost(double traversal_cost, double intersection_cost) const;
  // Returns the smallest axis aligned bounding box of *this
  Voxel boundingBox() const;

 private:
  std::unique_ptr<Node> root;
};

// Assumes `k_t` > 0 and `k_i` > 0
Tree buildKdTree(const Vector<SceneTriangle*>& scene_triangles, double k_t,
                 double k_i);
Vector<Triangle> extractTriangles(
    const Vector<SceneTriangle*>& scene_triangles);

// Number of triangles left, on and right of the sweeping plane
// in the kd tree building algorithm
struct TriangleCounts {
  // number of triagles having an overlap with non-zero area
  // with V_l \ p and V_r \p
  int left = 0;
  // ... with p
  int plane = 0;
  // ... with V_r \ p
  int right = 0;
};
// Calculates the surface area heuristic cost
// for splits [(n_left + n_plane), (n_right)]
// and [(n_left), (n_right + n_plane)]

// returns the cost of the better split and
// the side to which n_plane was assigned
// in this split
std::pair<double, bool> surfaceAreaHeuristic(double l_area, double r_area,
                                             TriangleCounts counts,
                                             double traversal_cost,
                                             double intersection_cost);
// Returns the areas of subvoxels formed by `plane` from `voxel`
// as a fraction of the area of `voxel`
std::pair<double, double> relativeSubvoxelAreas(const Voxel& voxel,
                                                const AxisPlane& plane);
#endif
