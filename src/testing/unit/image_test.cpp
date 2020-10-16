#include "../../image.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fstream>

using ::testing::ElementsAreArray;
namespace {
// Tests that the pixel values are stored at the correct
// positions of the file.
TEST(TestImage, Simple) {
  Image image(2, 3);
  image.setColor(1, 2, Vec3(1.0, 1.0, 0.0));
  image.setColor(0, 0, Vec3(1.0, 0.0, 1.0));
  image.savePPM(".test.ppm");
  std::ifstream fin(".test.ppm");
  std::vector<unsigned char> file;

  while (fin.good()) {
    unsigned char x = fin.get();
    if (fin.good()) {
      file.push_back(x);
    }
  }
  std::vector<unsigned char> correct_file = {
      'P',  '6', '\n', '2', ' ', '3', '\n', '2', '5', '5',
      '\n', 255, 0,    255, 0,   0,   0,    0,   0,   0,
      0,    0,   0,    0,   0,   0,   255,  255, 0};

  EXPECT_THAT(file, ElementsAreArray(correct_file));
}
// tests that values are correctly transformed
// from linear color space to Srgb
TEST(TestImage, LinearToSrgb) {
  Image image(1, 1);
  image.setColor(0, 0, Vec3(0.0, 0.5, 1.0));
  image.savePPM(".test.ppm");
  std::ifstream fin(".test.ppm");
  std::vector<unsigned char> file;

  while (fin.good()) {
    unsigned char x = fin.get();
    if (fin.good()) {
      file.push_back(x);
    }
  }
  // 0.5 * 255 in linear color space corresponds
  // to about 188 in Srgb color space
  std::vector<unsigned char> correct_file = {
      'P',  '6', '\n', '1', ' ', '1', '\n', '2', '5', '5',
      '\n', 0, 188,  255};

  EXPECT_THAT(file, ElementsAreArray(correct_file));
}
TEST(TestImage, ScaleMaxTo) {
  Image image(2, 3);
  image.setColor(1, 2, Vec3(10, 5, 0.0));
  image.setColor(0, 0, Vec3(10.0, 0.0, 10.0));
  image.scaleMaxTo(1.0);
  image.savePPM(".test.ppm");
  std::ifstream fin(".test.ppm");
  std::vector<unsigned char> file;

  while (fin.good()) {
    unsigned char x = fin.get();
    if (fin.good()) {
      file.push_back(x);
    }
  }
  std::vector<unsigned char> correct_file = {
      'P',  '6', '\n', '2', ' ', '3', '\n', '2', '5', '5',
      '\n', 255, 0,    255, 0,   0,   0,    0,   0,   0,
      0,    0,   0,    0,   0,   0,   255,  188, 0};

  EXPECT_THAT(file, ElementsAreArray(correct_file));
}
TEST(TestImage, DistanceTo) {
  Image image(2, 3);
  Image image2(2, 3);
  image.setColor(1, 2, Vec3(10, 5, 0.0));
  image2.setColor(1, 2, Vec3(10, 5, 0.0));
  EXPECT_NEAR(image.distanceTo(image2), 0.0, EPS);
  image.setColor(1, 0, Vec3(1.0, 0.0, 0.0));
  EXPECT_NEAR(image.distanceTo(image2), 1.0, EPS);
  image2.setColor(1, 1, Vec3(0.0, 0.0, 2.0));
  EXPECT_NEAR(image.distanceTo(image2), std::sqrt(1.0 + 4.0), EPS);
}
}  // namespace
