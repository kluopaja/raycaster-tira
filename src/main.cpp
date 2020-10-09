#include <iostream>

#include "raycaster.h"
#include "image.h"
void renderBallScene() {
  Camera camera(Vec3(0, 2, 5.0), Vec3(0.0, -0.2, -1.0),
                Vec3(0.0, 1.0, -0.2), kPi / 2.0, kPi / 2.0);
  Scene scene(camera);
  scene.setEnvironmentLightColor(Vec3(0.2, 0.2, 0.2));
  scene.addModelFromFile("../models/large_white_light.obj", Vec3(20, 10.0, 5.0), kSmooth);
  scene.addPointLight(Vec3(-10.0, 20.0, -10.0), Vec3(3.0));
  //scene.setSamplingScheme(kUniformSphere);
  scene.addModelFromFile("../models/scene.obj", Vec3(0.0), kSmooth);
  Image result = scene.render(500, 500, 10, 4, 4);
  result.scaleMaxTo(1.0);
  result.savePPM("test.ppm");
  result.truncateToFraction(0.99);
  result.scaleMaxTo(1.0);
  result.savePPM("test_truncated.ppm");
  result.toLog2(0.000001);
  result.scaleMaxTo(1.0);
  result.savePPM("test_log.ppm");
}
int main() {
  // this image can be used to check if your image viewer can
  // display the PPM files correctly
  generateGammaTestImage(1000, 1000).savePPM("gamma_test.ppm");
  renderBallScene();
}
