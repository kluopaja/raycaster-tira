#include "geometry.h"
#include <random>
Vec3 randomVec3(double lo, double hi, std::mt19937& random_engine);
double randomLogUniformReal(double lo_log, double hi_log,
                            std::mt19937& random_engine);
Vec2 randomBaryCoords(std::mt19937& random_engine);
std::vector<Triangle*> randomTriangleVector(double lo, double hi, int n,
                                            std::mt19937 random_engine);
