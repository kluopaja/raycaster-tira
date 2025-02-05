#include "kd_tree.h"

#include <limits>
#include <memory>
#include <tuple>

#include "raycaster.h"
#include "utils.h"
#include "vector.h"
namespace {

class Event {
 public:
  int type;  // 0 == +, 1 == |, 2 == -
  double pos;
};
// class for storing the triangle during the build process
struct BuildTriangle {
  ClipTriangle clip_triangle;
  SceneTriangle* scene_triangle;
};
Vector<BuildTriangle> createBuildTriangles(
    const Vector<SceneTriangle*>& scene_triangles) {
  Vector<BuildTriangle> build_triangles;
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    build_triangles.pushBack(
        {ClipTriangle(scene_triangles[i]->triangle), scene_triangles[i]});
  }
  return build_triangles;
}
Vector<SceneTriangle*> extractSceneTriangles(
    const Vector<BuildTriangle>& build_triangles) {
  Vector<SceneTriangle*> scene_triangles;
  for (size_t i = 0; i < build_triangles.size(); ++i) {
    scene_triangles.pushBack(build_triangles[i].scene_triangle);
  }
  return scene_triangles;
}
bool operator<(const Event& a, const Event& b) { return a.pos < b.pos; }
struct SplitPlane {
  AxisPlane plane;
  bool side;
  double cost;
};
struct SplitResult {
  Vector<BuildTriangle> left;
  Vector<BuildTriangle> right;
};
class TreeBuilder {
 public:
  TreeBuilder(double traversal_cost, double intersection_cost);
  Tree build(const Vector<SceneTriangle*>& scene_triangles) const;

 private:
  double traversal_cost;
  double intersection_cost;

  std::unique_ptr<Node> recursiveBuild(Vector<BuildTriangle>& build_triangles,
                                       const Voxel& voxel) const;
  SplitPlane findPlane(const Vector<BuildTriangle>& clip_triangles,
                       const Voxel& v) const;

  Vector<Event> createEventList(const Vector<BuildTriangle>& build_triangles,
                                int dimension) const;

  SplitResult splitTriangles(const Vector<BuildTriangle>& build_triangles,
                             const SplitPlane& plane) const;
};
TreeBuilder::TreeBuilder(double traversal_cost, double intersection_cost)
    : traversal_cost(traversal_cost), intersection_cost(intersection_cost) {
  assert(traversal_cost > EPS && intersection_cost > EPS);
}

Tree TreeBuilder::build(const Vector<SceneTriangle*>& scene_triangles) const {
  Voxel voxel = boundingBox(extractTriangles(scene_triangles));
  Vector<BuildTriangle> build_triangles = createBuildTriangles(scene_triangles);
  std::unique_ptr<Node> root = recursiveBuild(build_triangles, voxel);
  return Tree(std::move(root));
}
// implements the O(n log^2 n) algorithm
std::unique_ptr<Node> TreeBuilder::recursiveBuild(
    Vector<BuildTriangle>& build_triangles, const Voxel& voxel) const {
  if (build_triangles.size() == 0) {
    return std::make_unique<Node>(Vector<Triangle>(), Vector<SceneTriangle*>(),
                                  voxel);
  }
  // the surface area heuristic is not well defined if the voxel has
  // zero area so we will just terminate
  if (voxel.area() < EPS) {
    // make a leaf
    Vector<SceneTriangle*> scene_triangles =
        extractSceneTriangles(build_triangles);
    return std::make_unique<Node>(extractTriangles(scene_triangles),
                                  scene_triangles, voxel);
  }
  SplitPlane split_plane = findPlane(build_triangles, voxel);
  // termination criterion
  if (split_plane.cost > intersection_cost * (double)build_triangles.size()) {
    // make a leaf
    Vector<SceneTriangle*> scene_triangles =
        extractSceneTriangles(build_triangles);
    return std::make_unique<Node>(extractTriangles(scene_triangles),
                                  scene_triangles, voxel);
  }
  SplitResult split = splitTriangles(build_triangles, split_plane);
  build_triangles.clear();
  Voxel left_voxel(voxel);
  Voxel right_voxel(voxel);
  left_voxel.clip(split_plane.plane, 0);
  right_voxel.clip(split_plane.plane, 1);
  std::unique_ptr<Node> left = recursiveBuild(split.left, left_voxel);
  std::unique_ptr<Node> right = recursiveBuild(split.right, right_voxel);
  return std::make_unique<Node>(std::move(left), std::move(right), voxel,
                                split_plane.plane);
}
SplitPlane TreeBuilder::findPlane(const Vector<BuildTriangle>& build_triangles,
                                  const Voxel& voxel) const {
  SplitPlane best_split = {{}, 0, std::numeric_limits<double>::infinity()};
  for (int i = 0; i < 3; ++i) {
    Vector<Event> event_list = createEventList(build_triangles, i);
    TriangleCounts counts;
    counts.left = 0;
    counts.plane = 0;
    counts.right = build_triangles.size();
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
      // _\ .
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
      int n_new_on_plane = 0;
      while (j < event_list.size() &&
             event_list[j].pos <= current_plane.pos + EPS) {
        if (event_list[j].type == 0) {
          ++n_start_plane_intersection;
        }
        if (event_list[j].type == 1) {
          ++n_new_on_plane;
        }
        if (event_list[j].type == 2) {
          ++n_end_plane_intersection;
        }
        ++j;
      }
      // updating intersection counts
      counts.right -= n_end_plane_intersection;
      counts.right -= n_new_on_plane;
      counts.plane += n_new_on_plane;
      // Plane check
      double l_area, r_area;
      std::tie(l_area, r_area) = relativeSubvoxelAreas(voxel, current_plane);
      double cost;
      bool side;
      std::tie(cost, side) = surfaceAreaHeuristic(
          l_area, r_area, counts, traversal_cost, intersection_cost);
      if (cost < best_split.cost) {
        best_split = {current_plane, side, cost};
      }
      // After
      counts.left += n_start_plane_intersection;
      counts.plane -= n_new_on_plane;
      counts.left += n_new_on_plane;
    }
  }
  return best_split;
}
Vector<Event> TreeBuilder::createEventList(
    const Vector<BuildTriangle>& build_triangles, int dimension) const {
  Vector<Event> event_list;
  for (size_t i = 0; i < build_triangles.size(); ++i) {
    if (build_triangles[i].clip_triangle.isAxisAligned(dimension)) {
      event_list.pushBack({1, build_triangles[i].clip_triangle.min(dimension)});
    } else {
      event_list.pushBack({0, build_triangles[i].clip_triangle.min(dimension)});
      event_list.pushBack({2, build_triangles[i].clip_triangle.max(dimension)});
    }
  }
  quickSort(event_list.begin(), event_list.end());
  return event_list;
}
SplitResult TreeBuilder::splitTriangles(
    const Vector<BuildTriangle>& build_triangles,
    const SplitPlane& plane) const {
  // Split build_triangles according to the plane
  Vector<BuildTriangle> left_build_triangles;
  Vector<BuildTriangle> right_build_triangles;
  for (size_t i = 0; i < build_triangles.size(); ++i) {
    bool to_left, to_right;
    std::tie(to_left, to_right) =
        build_triangles[i].clip_triangle.overlapsSides(plane.plane, plane.side);
    if (!to_left) {
      right_build_triangles.pushBack(build_triangles[i]);
    } else if (!to_right) {
      left_build_triangles.pushBack(build_triangles[i]);
    } else {
      BuildTriangle left = build_triangles[i];
      BuildTriangle right = build_triangles[i];
      left.clip_triangle.clip(plane.plane, 0);
      right.clip_triangle.clip(plane.plane, 1);
      left_build_triangles.pushBack(left);
      right_build_triangles.pushBack(right);
    }
  }
  return {left_build_triangles, right_build_triangles};
}

}  // namespace

Node::Node(std::unique_ptr<Node> left, std::unique_ptr<Node> right,
           const Voxel& voxel, const AxisPlane& plane)
    : left(std::move(left)),
      right(std::move(right)),
      voxel(voxel),
      plane(plane) {}
Node::Node(const Vector<Triangle>& triangles,
           const Vector<SceneTriangle*>& scene_triangles, const Voxel& voxel)
    : triangles(triangles), scene_triangles(scene_triangles), voxel(voxel) {}

ScenePoint Node::getClosestRayIntersection(const Ray& r) const {
  if (!voxel.intersects(r)) {
    return {nullptr, {}};
  }
  if (isLeaf()) {
    RayTriangleIntersection res = firstRayTriangleIntersection(triangles, r);
    // Check that the intersection is actually within the voxel.
    // Some of the triangles extend outside the voxels
    // but the tree traversal assumes that the intersections are within
    // the voxel (the closer voxel is checked first for the intersections)
    if (res.index != triangles.size() &&
        voxel.isInside(triangles[res.index].pointFromBary(res.bary_coords))) {
      return {scene_triangles[res.index], res.bary_coords};
    }
    return {nullptr, {}};
  }
  // first the closer leaf
  ScenePoint result;
  if (r.getDirection()[plane.axis] >= 0) {
    result = left->getClosestRayIntersection(r);
    if (result.scene_triangle != nullptr) {
      return result;
    }
    return right->getClosestRayIntersection(r);
  }
  result = right->getClosestRayIntersection(r);
  if (result.scene_triangle != nullptr) {
    return result;
  }
  return left->getClosestRayIntersection(r);
}
std::size_t Node::numLeaves() const {
  if (isLeaf()) return 1;
  return left->numLeaves() + right->numLeaves();
}
std::size_t Node::numNonEmptyLeaves() const {
  if (isLeaf()) {
    if (triangles.size() > 0) return 1;
    return 0;
  }
  return left->numNonEmptyLeaves() + right->numNonEmptyLeaves();
}
std::size_t Node::totalTriangles() const {
  if (isLeaf()) {
    return triangles.size();
  }
  return left->totalTriangles() + right->totalTriangles();
}
double Node::cost(double traversal_cost, double intersection_cost) const {
  double cost_sum = 0;
  if (voxel.area() < EPS) return 0;
  // I assume that they didn't include the leaves to this
  if (!isLeaf()) {
    cost_sum += voxel.area() * traversal_cost;
  }
  if (isLeaf()) {
    cost_sum += voxel.area() * triangles.size() * intersection_cost;
    return cost_sum / voxel.area();
  }
  cost_sum +=
      left->cost(traversal_cost, intersection_cost) * left->voxel.area();
  cost_sum +=
      right->cost(traversal_cost, intersection_cost) * right->voxel.area();
  return cost_sum / voxel.area();
};
Voxel Node::boundingBox() const { return voxel; }
std::size_t Tree::numLeaves() const { return root->numLeaves(); }
std::size_t Tree::numNonEmptyLeaves() const {
  return root->numNonEmptyLeaves();
}
double Tree::averageTrianglesInLeaf() const {
  return (double)root->totalTriangles() / root->numNonEmptyLeaves();
}
double Tree::cost(double traversal_cost, double intersection_cost) const {
  return root->cost(traversal_cost, intersection_cost);
};
Voxel Tree::boundingBox() const { return root->boundingBox(); }
ScenePoint Tree::getClosestRayIntersection(const Ray& r) const {
  return root->getClosestRayIntersection(r);
};
// for now simply uses the getClosestRayIntersection
// could possibly be sped up implementing a seconds recursive
// intersection test separately for this
bool Tree::trianglesIntersectSegment(const Vec3& a, const Vec3& b) const {
  Ray r(a, b - a);
  ScenePoint sp = root->getClosestRayIntersection(r);
  // no intersection
  if (sp.scene_triangle == nullptr) {
    return false;
  }
  Vec3 intersection_point =
      sp.scene_triangle->triangle.pointFromBary(sp.bary_coords);
  if (pointOnSegment(intersection_point, a, b)) {
    return true;
  }
  return false;
}
Vector<Triangle> extractTriangles(
    const Vector<SceneTriangle*>& scene_triangles) {
  Vector<Triangle> triangles;
  for (size_t i = 0; i < scene_triangles.size(); ++i) {
    triangles.pushBack(scene_triangles[i]->triangle);
  }
  return triangles;
}
// Favors splits that cut off empty space
//
// Corresponds to \lambda in the paper
inline double cutOffBonus(double l_area, double r_area, TriangleCounts counts) {
  // Check that some empty space is really being cut off
  if (l_area < 1.0 - EPS && counts.right == 0) {
    return 0.8;
  }
  if (r_area < 1.0 - EPS && counts.left == 0) {
    return 0.8;
  }
  return 1.0;
}
std::pair<double, bool> surfaceAreaHeuristic(double l_area, double r_area,
                                             TriangleCounts counts,
                                             double traversal_cost,
                                             double intersection_cost) {
  std::pair<double, bool> left;
  // try inserting triangles on plane to the left subtree
  double cut_off_bonus = cutOffBonus(l_area, r_area, counts);
  left.first = traversal_cost +
               intersection_cost * (l_area * (counts.left + counts.plane) +
                                    r_area * counts.right);

  left.first *= cut_off_bonus;
  left.second = false;
  std::pair<double, bool> right;
  right.first = traversal_cost +
                intersection_cost * (l_area * counts.left +
                                     r_area * (counts.right + counts.plane));
  right.first *= cut_off_bonus;
  right.second = true;
  return std::min(left, right);
}
// returns SA(V_l)/SA(V) and SA(V_r)/SA(V)
std::pair<double, double> relativeSubvoxelAreas(const Voxel& voxel,
                                                const AxisPlane& plane) {
  // axes for the voxel face parallel to plane
  int base1 = (plane.axis + 1) % 3;
  int base2 = (plane.axis + 2) % 3;
  double base1_length = voxel.hi[base1] - voxel.lo[base1];
  double base2_length = voxel.hi[base2] - voxel.lo[base2];
  double side_length = voxel.hi[plane.axis] - voxel.lo[plane.axis];

  double base_area = base1_length * base2_length;
  double base_perimeter = 2 * (base1_length + base2_length);
  double side_area = base_perimeter * side_length;
  double left_length = (plane.pos - voxel.lo[plane.axis]);

  double sa_left = base_perimeter * left_length + 2 * base_area;
  double sa_right =
      base_perimeter * (side_length - left_length) + 2 * base_area;
  double sa_total = side_area + 2 * base_area;
  assert(sa_total > EPS);
  return {sa_left / sa_total, sa_right / sa_total};
}
Tree buildKdTree(const Vector<SceneTriangle*>& triangles, double k_t,
                 double k_i) {
  assert(k_t > EPS && k_i > EPS);
  TreeBuilder builder(k_t, k_i);
  return builder.build(triangles);
}
