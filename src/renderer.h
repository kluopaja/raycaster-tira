#ifndef RAYCASTER_RENDERER_H
#define RAYCASTER_RENDERER_H
#include <string>

#include "input.h"
void executeRecipe(const Recipe& recipe);
void executeInputFile(const std::string& file_name, bool verbose);
#endif
