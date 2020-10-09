#ifndef RAYCASTER_MODEL_LOADER_H
#define RAYCASTER_MODEL_LOADER_H
#include <string>

#include "raycaster.h"

bool loadModel(const std::string& file, Model& out_model,
               NormalType normal_type);
#endif
