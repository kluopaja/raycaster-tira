#ifndef RAYCASTER_IMAGE_H
#define RAYCASTER_IMAGE_H
#include <cmath>
#include <string>

#include "geometry.h"
class Image {
 public:
  Image(int x_resolution, int y_resolution);
  void setColor(int x_pos, int y_pos, const Vec3& color);
  Vec3 getColor(int x_pos, int y_pos);
  // scales values to [0, 1]
  void scaleMaxTo(double new_max);
  // applies value = log_2((value + c)/c)
  // for debugging.
  void toLog2(double c);
  // Let f be smallest value such that c of all pixel values are at most f
  // then applies value = max(value, f)
  void truncateToFraction(double c);
  // saves image as PPM with sRGB encoded values
  bool savePPM(const std::string& file);
  // Calculates the Euclidean distance between two images
  double distanceTo(const Image& other);
 private:
  double linearToSrgb(double val);
  double maxColorValue();
  int bufferPos(int x, int y);
  char floatToByte(double val);
  int x_resolution;
  int y_resolution;
  std::vector<Vec3> buffer;
};
// can be used to check if the used ppm viewer displays the image correctly
// the whole image should look gray when looked at
// display native resolution
// (if separate black and white pixels are seen on the left side, then
// the monitor resolution is too small or the monitor is being looked
// from a too short distance)
Image generateGammaTestImage(int x_size, int y_size);
inline void Image::setColor(int x, int y, const Vec3& color) {
  assert(x >= 0 && x < x_resolution);
  assert(y >= 0 && y < y_resolution);
  assert(color[0] > -EPS && color[1] > -EPS && color[2] > -EPS);
  buffer[(size_t)bufferPos(x, y)] = color;
}
inline Vec3 Image::getColor(int x, int y) {
  return buffer[(size_t)bufferPos(x, y)];
}
inline int Image::bufferPos(int x, int y) { return y * x_resolution + x; }
inline char Image::floatToByte(double val) { return std::round(val * 255); }
#endif
