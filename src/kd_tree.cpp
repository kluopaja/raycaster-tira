#include "kd_tree.h"

#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

Node::Node(Node* left, Node* right, const Voxel& voxel, const AxisPlane& plane)
    : left(left), right(right), voxel(voxel), plane(plane) {}
Node::Node(const Voxel& voxel, const std::vector<Triangle*>& triangles)
    : voxel(voxel), triangles(triangles) {}
Node::~Node() {
  delete left;
  delete right;
}
TrianglePoint Node::getClosestRayIntersection(const Ray& r) const {
  if (!voxel.intersects(r)) {
    return {nullptr, {}};
  }
  if (isLeaf()) {
    return getClosestRayLeafIntersection(r);
  }
  // first the closer leaf
  TrianglePoint result;
  if (r.direction[plane.axis] >= 0) {
    result = left->getClosestRayIntersection(r);
    if (result.triangle != nullptr) {
      return result;
    }
    return right->getClosestRayIntersection(r);
  }
  result = right->getClosestRayIntersection(r);
  if (result.triangle != nullptr) {
    return result;
  }
  return left->getClosestRayIntersection(r);
}
Tree::Tree(Tree&& a) noexcept : root(root) { a.root = nullptr; }
TrianglePoint Tree::getClosestRayIntersection(const Ray& r) const {
  return root->getClosestRayIntersection(r);
};
bool operator<(const Event& a, const Event& b) { return a.pos < b.pos; }
TreeBuilder::TreeBuilder(double traversal_cost, double intersection_cost)
    : traversal_cost(traversal_cost), intersection_cost(intersection_cost) {}
Tree TreeBuilder::build(const std::vector<Triangle*>& triangles) const {
  Voxel voxel = boundingBox(triangles);
  std::vector<ClipTriangle> clip_triangles = createClipTriangles(triangles);
  Node* root = recursiveBuild(clip_triangles, voxel);
  return Tree(root);
}
std::vector<ClipTriangle> createClipTriangles(
    const std::vector<Triangle*>& triangles) const {
  std::vector<ClipTriangle> clip_triangles;
  for (size_t i = 0; i < triangles.size(); ++i) {
    clip_triangles.emplace_back(triangles[i]);
  }
  return clip_triangles;
}
std::vector<Triangle*> extractTriangles(
    const std::vector<ClipTriangle>& clip_triangles) const {
  std::vector<Triangle*> triangles;
  for (size_t i = 0; i < clip_triangles.size(); ++i) {
    triangles.push_back(clip_triangles[i].triangle);
  }
  return triangles;
}
// implements the O(n log^2 n) algorithm
Node* TreeBuilder::recursiveBuild(std::vector<ClipTriangle>& clip_triangles,
                                  const Voxel& voxel) const {
  SplitPlane split_plane = findPlane(clip_triangles, voxel);
  // termination criterion
  if (split_plane.cost > intersection_cost * (double)clip_triangles.size()) {
    // make a leaf
    return new Node(voxel, extractTriangles(clip_triangles));
  }
  SplitResult split = splitTriangles(clip_triangles, split_plane);
  clip_triangles.clear();
  Voxel left_voxel(voxel);
  Voxel right_voxel(voxel);
  left_voxel.clip(split_plane.plane, 0);
  right_voxel.clip(split_plane.plane, 1);
  Node* left = recursiveBuild(split.left, left_voxel);
  Node* right = recursiveBuild(split.right, right_voxel);
  Node* node = new Node(left, right, voxel, split_plane.plane);
  return node;
}
SplitPlane TreeBuilder::findPlane(
    const std::vector<ClipTriangle>& clip_triangles, const Voxel& voxel) const {
  SplitPlane best_split = {{}, 0, std::numeric_limits<double>::infinity()};
  for (int i = 0; i < 3; ++i) {
    std::vector<Event> event_list = createEventList(clip_triangles, i);
    // number of triagles having an overlap with non-zero area
    // with V_l \ p and V_r \p
    int n_left = 0;
    // ... with p
    int n_plane = 0;
    // ... with V_r \ p
    int n_right = clip_triangles.size();
    AxisPlane current_plane = {i, 0.0};
    for (size_t j = 0; j < event_list.size();) {
      current_plane.pos = event_list[j].pos;
      //  :/
      //  /_
      //  :
      // Before plane check:
      //     These are still counted as part of the right side
      // After:
      //     Added to left side
      int n_start_plane_intersection = 0;
      // \:
      // _\
            //  :
      // Before plane check:
      //     Removed from the right side
      int n_end_plane_intersection = 0;
      //  :
      //  |
      //  :
      // Before plane check:
      //     Removed from the right side
      //     Added to the plane
      // After:
      //     Added to left side
      int n_on_plane = 0;
      while (j < event_list.size() &&
             event_list[j].pos <= current_plane.pos + EPS) {
        if (event_list[j].type == 0) {
          ++n_end_plane_intersection;
        }
        if (event_list[j].type == 1) {
          ++n_on_plane;
        }
        if (event_list[j].type == 2) {
          ++n_start_plane_intersection;
        }
        ++j;
      }
      n_right -= n_end_plane_intersection;
      n_right -= n_on_plane;
      // Plane check
      double l_area, r_area;
      std::tie(l_area, r_area) = relative_subvoxel_areas(voxel, current_plane);
      double cost;
      bool side;
      std::tie(cost, side) =
          surface_area_heuristic(l_area, r_area, n_left, n_plane, n_right);
      if (cost < best_split.cost) {
        best_split = {current_plane, side, cost};
      }
      // After
      n_left += n_start_plane_intersection;
      n_left += n_on_plane;
    }
  }
  return best_split;
}
std::vector<Event> TreeBuilder::createEventList(
    const std::vector<ClipTriangle>& clip_triangles, int dimension) const {
  std::vector<Event> event_list;
  for (size_t i = 0; i < clip_triangles.size(); ++i) {
    if (clip_triangles[i].isAxisAligned(dimension)) {
      event_list.push_back({1, clip_triangles[i].min(dimension)});
    } else {
      event_list.push_back({0, clip_triangles[i].min(dimension)});
      event_list.push_back({2, clip_triangles[i].max(dimension)});
    }
  }
  std::sort(event_list.begin(), event_list.end());
  return event_list;
}
SplitResult TreeBuilder::splitTriangles(
    const std::vector<ClipTriangle>& clip_triangles,
    const SplitPlane& plane) const {
  // Split clip_triangles according to the plane
  std::vector<ClipTriangle> left_clip_triangles;
  std::vector<ClipTriangle> right_clip_triangles;
  for (int i = 0; i < clip_triangles.size(); ++i) {
    bool to_left, to_right;
    std::tie(to_left, to_right) =
        clip_triangles[i].overlapsSides(plane.plane, plane.side);
    if (!to_left) {
      right_clip_triangles.push_back(clip_triangles[i]);
    } else if (!to_right) {
      left_clip_triangles.push_back(clip_triangles[i]);
    } else {
      ClipTriangle left = clip_triangles[i];
      ClipTriangle right = clip_triangles[i];
      left.clip(plane.plane, 0);
      right.clip(plane.plane, 1);
      left_clip_triangles.push_back(left);
      right_clip_triangles.push_back(right);
    }
  }
  return {left_clip_triangles, right_clip_triangles};
}

std::pair<double, bool> surfaceAreaHeuristic(
    double l_area, double r_area, int n_left, int n_plane, int n_right) const {
  std::pair<double, bool> left;
  // try inserting triangles on plane to the left subtree
  left.first =
      traversal_cost +
      intersection_cost * (l_area * (n_left + n_plane) + r_area * n_right);
  left.second = false;
  std::pair<double, bool> right;
  right.first =
      traversal_cost +
      intersection_cost * (l_area * n_left + r_area * (n_right + n_plane));
  right.second = true;
  return std::min(left, right);
}
// returns SA(V_l)/SA(V) and SA(V_l)/SA(V)
std::pair<double, double> relativeSubvoxelAreas(
    const Voxel& voxel, const AxisPlane& plane) const {
  // axes for the voxel face parallel to plane
  int base1 = (plane.axis + 1) % 3;
  int base2 = (plane.axis + 2) % 3;
  double base1_length = voxel.hi[base1] - voxel.lo[base1];
  double base2_length = voxel.hi[base2] - voxel.lo[base2];
  double side_length = voxel.hi[plane.axis] - voxel.lo[plane.axis];

  double base_area = base1_length * base2_length;
  double base_perimeter = base1_length + base2_length;
  double side_area = base_perimeter * side_length;
  double left_length = (plane.pos - voxel.lo[plane.axis]);

  double sa_left = base_perimeter * left_length + base_area;
  double sa_right = base_perimeter * (side_length - left_length) + base_area;
  double sa_total = side_area + base_area;
  return {sa_left / sa_total, sa_right / sa_total};
}
Tree buildKdTree(std::vector<Triangle*> triangles, double k_t, double k_i) {
  TreeBuilder builder(k_t, k_i);
  return builder.build(triangles);
}
