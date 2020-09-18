#include "geometry.h"
#ifndef RAYCASTER_KD_TREE_H
#define RAYCASTER_KD_TREE_H

class Node {
 public:
  Node(Node* left, Node* right, const Voxel& voxel, const AxisPlane& plane);
  Node(const Voxel& voxel, const std::vector<Triangle*>& triangles);
  ~Node();
  TrianglePoint getClosestRayIntersection(const Ray& r) const;

 private:
  // the side with smaller coordinates
  // for example if the plane axis = 0,
  // then maximum x coordinate of left.voxel is at most the
  // minimum x coordinate of right.voxel
  Node* left;
  Node* right;
  std::vector<Triangle*> triangles;
  Voxel voxel;
  AxisPlane plane;
  bool isLeaf() const;
};
inline bool Node::isLeaf() const {
  return (left == nullptr && right == nullptr);
}
class Tree {
 public:
  TrianglePoint getClosestRayIntersection(const Ray& r) const;
  Tree(Node* root) : root(root) {}
  // move constructor
  Tree(Tree&& a) noexcept;
  ~Tree() { delete root; }

 private:
  Node* root;
};

Tree buildKdTree(const std::vector<Triangle*>& triangles, double k_t,
                 double k_i);

std::vector<ClipTriangle> createClipTriangles(
    const std::vector<Triangle*>& triangles);
std::vector<Triangle*> extractTriangles(
    const std::vector<ClipTriangle>& clip_triangles);
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
