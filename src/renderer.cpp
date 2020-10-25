#include "renderer.h"

#include <string>

#include "input.h"
#include "raycaster.h"

void executeRecipe(const Recipe& recipe) {
  Camera camera(recipe.camera.pos, recipe.camera.front, recipe.camera.up,
                recipe.camera.x_fov, recipe.camera.y_fov);
  Scene scene(camera);
  for (auto& model : recipe.models) {
    scene.addModelFromFile(model.file, model.pos, model.normal);
  }
  for (auto& point_light : recipe.point_lights) {
    scene.addPointLight(point_light.pos, point_light.color);
  }
  scene.setEnvironmentLightColor(recipe.environment_light.color);
  if (recipe.environment_light.type == EnvironmentLightType::kDirected) {
    scene.setEnvironmentLightDirected(recipe.environment_light.direction,
                                      recipe.environment_light.exp);
  }
  Image image = scene.render(recipe.render.width, recipe.render.height,
                             recipe.render.pixel_rays, recipe.render.depth,
                             recipe.render.branching);
  for (auto& config : recipe.images) {
    Image image_copy = image;
    image.truncateToFraction(config.truncate);
    image.scaleMaxTo(config.scale_max);
    image.savePPM(config.file);
  }
}
void executeInputFile(const std::string& file_name, bool verbose) {
  Recipe recipe = loadRecipe(file_name, verbose);
  executeRecipe(recipe);
}
