#include "image.h"

#include <fstream>

#include "utils.h"

Image::Image(int x_resolution, int y_resolution)
    : x_resolution(x_resolution), y_resolution(y_resolution) {
  assert(x_resolution > 0);
  assert(y_resolution > 0);
  buffer.resize((size_t)x_resolution * (size_t)y_resolution);
}
void Image::scaleMaxTo(double new_max) {
  assert(new_max > EPS);
  double old_max = maxColorValue();
  if (old_max < EPS) {
    return;
  }
  double scale_factor = new_max / old_max;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      size_t pos = (size_t)bufferPos(x, y);
      buffer[pos] = buffer[pos] * scale_factor;
    }
  }
}
void Image::truncateToFraction(double c) {
  assert(c > 0);
  assert(c < 1 + EPS);
  std::vector<double> values;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      for (int i = 0; i < 3; ++i) {
        values.push_back(buffer[(size_t)bufferPos(x, y)][i]);
      }
    }
  }
  quickSort(values.begin(), values.end());
  assert(values.size() > 0);
  std::size_t fraction_pos =
      std::max((std::size_t)1, (std::size_t)(c * values.size() + EPS)) - 1.0;
  double limit = values[fraction_pos];
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      for (int i = 0; i < 3; ++i) {
        size_t pos = (size_t)bufferPos(x, y);
        buffer[pos][i] = std::min(buffer[pos][i], limit);
      }
    }
  }
}
void Image::toLog2(double c) {
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      size_t pos = (size_t)bufferPos(x, y);
      // r g b
      for (int i = 0; i < 3; ++i) {
        buffer[pos][i] = std::log2(buffer[pos][i] + c) - std::log2(c);
      }
    }
  }
}
bool Image::savePPM(const std::string& file) {
  std::ofstream fout(file, std::ios::binary);
  if (!fout.is_open()) {
    return false;
  }
  // binary mode
  fout << "P6" << std::endl;
  fout << x_resolution << " " << y_resolution << std::endl;
  fout << "255" << std::endl;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      // r g b
      for (int i = 0; i < 3; ++i) {
        fout << floatToByte(linearToSrgb(buffer[(size_t)bufferPos(x, y)][i]));
      }
    }
  }
  fout.close();
  return true;
}
double Image::distanceTo(const Image& other) {
  assert(x_resolution == other.x_resolution);
  assert(y_resolution == other.y_resolution);
  double sum_of_squares = 0;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      // r g b
      for (int i = 0; i < 3; ++i) {
        std::size_t pos = (size_t)bufferPos(x, y);
        sum_of_squares += std::pow(buffer[pos][i] - other.buffer[pos][i], 2.0);
      }
    }
  }
  return std::sqrt(sum_of_squares);
}
// Converts linear intensities to sRGB values
// ("gamma" encoding)
// Human eye perceives light in a logarithmic scale.
// Therefore encoding the light intensities without
// any transformation would leave little (perceived)
// resolution to dim intensities.
//
// Note that this is not a logarithmic tranformation but
// apparently still works well.
//
// see http://color.org/chardata/rgb/sRGB.pdf
double Image::linearToSrgb(double val) {
  if (val <= 0.0031308) {
    return 12.92 * val;
  }
  return 1.055 * std::pow(val, 1.0 / 2.4) - 0.055;
}
double Image::maxColorValue() {
  double max_value = 0;
  for (int y = 0; y < y_resolution; ++y) {
    for (int x = 0; x < x_resolution; ++x) {
      for (int i = 0; i < 3; ++i) {
        max_value = std::max(max_value, buffer[(size_t)bufferPos(x, y)][i]);
      }
    }
  }
  return max_value;
}
Image generateGammaTestImage(int x_size, int y_size) {
  assert(x_size > 20 && y_size > 20);
  Image image(x_size, y_size);
  for (int y = 0; y < y_size; ++y) {
    for (int x = 0; x < x_size / 2; ++x) {
      Vec3 color(0.0);
      if ((x + y) % 2 == 0) {
        color = Vec3(1.0);
      }
      image.setColor(x, y, color);
    }
    for (int x = x_size / 2; x < x_size; ++x) {
      Vec3 color(0.5);
      image.setColor(x, y, color);
    }
  }
  return image;
}
