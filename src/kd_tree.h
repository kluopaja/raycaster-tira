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
  TrianglePoint getClosestRayLeafIntersection(const Ray& r) const;
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
class Event {
 public:
  int type;  // 0 == +, 1 == |, 2 == -
  double pos;
};
bool operator<(const Event& a, const Event& b);
struct SplitPlane {
  AxisPlane plane;
  bool side;
  double cost;
};
struct SplitResult {
  std::vector<ClipTriangle> left;
  std::vector<ClipTriangle> right;
};
class TreeBuilder {
 public:
  TreeBuilder(double traversal_cost, double intersection_cost);
  Tree build(const std::vector<Triangle*>& triangles) const;

 private:
  double traversal_cost;
  double intersection_cost;

  Voxel boundingBox(const std::vector<Triangle*>& triangles) const;
  std::vector<ClipTriangle> createClipTriangles(
      const std::vector<Triangle*>& triangles) const;
  std::vector<Triangle*> extractTriangles(
      const std::vector<ClipTriangle>& clip_triangles) const;
  Node* recursiveBuild(std::vector<ClipTriangle>& clip_triangles,
                       const Voxel& voxel) const;
  SplitPlane findPlane(const std::vector<ClipTriangle>& clip_triangles,
                       const Voxel& v) const;

  std::vector<Event> createEventList(
      const std::vector<ClipTriangle>& clip_triangles, int dimension) const;

  SplitResult splitTriangles(const std::vector<ClipTriangle>& clip_triangles,
                             const SplitPlane& plane) const;
  // Calculates the surface area heuristic cost
  // for splits [(n_left + n_plane), (n_right)]
  // and [(n_left), (n_right + n_plane)]

  // returns the cost of the better split and
  // the side to which n_plane was assigned
  // in this split
  std::pair<double, bool> surface_area_heuristic(double l_area, double r_area,
                                                 int n_left, int n_plane,
                                                 int n_right) const;
  std::pair<double, double> relative_subvoxel_areas(
      const Voxel& voxel, const AxisPlane& plane) const;
};

Tree buildKdTree(const std::vector<Triangle*>& triangles, double k_t,
                 double k_i);

#endif
