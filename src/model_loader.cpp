#include "model_loader.h"

#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <assimp/Importer.hpp>

#include "geometry.h"
namespace {

class ModelLoader {
 public:
  bool load(const std::string& file, NormalType normal_type);
  Model model;

 private:
  void loadMaterials(const aiScene* ai_scene);
  void loadMaterial(const aiMaterial* ai_material);
  void loadNode(const aiNode* node, const aiScene* ai_scene);
  void loadMesh(const aiMesh* mesh);
};

bool ModelLoader::load(const std::string& file, NormalType normal_type) {
  Assimp::Importer importer;
  unsigned int normal_flag = aiProcess_GenNormals;
  if (normal_type == kSmooth) {
    normal_flag = aiProcess_GenSmoothNormals;
  }
  const aiScene* ai_scene = importer.ReadFile(
      file,
      aiProcess_Triangulate |  // every face will have <= 3 vertices
          normal_flag);
  if (!ai_scene) {
    std::cerr << "ERROR while loading the model: " << importer.GetErrorString()
              << std::endl;
    return false;
  }
  loadMaterials(ai_scene);
  loadNode(ai_scene->mRootNode, ai_scene);
  return true;
}
void ModelLoader::loadMaterials(const aiScene* ai_scene) {
  for (unsigned int i = 0; i < ai_scene->mNumMaterials; ++i) {
    loadMaterial(ai_scene->mMaterials[i]);
  }
}
void ModelLoader::loadMaterial(const aiMaterial* ai_material) {
  Material* new_material = new Material();
  aiColor3D diffuse(0.0, 0.0, 0.0);
  aiColor3D emitted(0.0, 0.0, 0.0);
  ai_material->Get(AI_MATKEY_COLOR_DIFFUSE, diffuse);
  ai_material->Get(AI_MATKEY_COLOR_AMBIENT, emitted);
  new_material->diffuse = Vec3(diffuse.r, diffuse.g, diffuse.b);
  new_material->emitted = Vec3(emitted.r, emitted.g, emitted.b);
  model.materials.push_back(new_material);
}
void ModelLoader::loadNode(const aiNode* node, const aiScene* ai_scene) {
  for (size_t i = 0; i < node->mNumMeshes; ++i) {
    aiMesh* mesh = ai_scene->mMeshes[node->mMeshes[i]];
    loadMesh(mesh);
  }
  for (size_t i = 0; i < node->mNumChildren; ++i) {
    loadNode(node->mChildren[i], ai_scene);
  }
}
void ModelLoader::loadMesh(const aiMesh* mesh) {
  for (size_t i = 0; i < mesh->mNumFaces; ++i) {
    SceneTriangle* scene_triangle = new SceneTriangle();
    aiFace face = mesh->mFaces[i];
    // ignore lines and points
    if (face.mNumIndices < 3) {
      continue;
    }
    assert(face.mNumIndices == 3);
    Vec3 vertices[3];
    for (int j = 0; j < 3; ++j) {
      vertices[j][0] = mesh->mVertices[face.mIndices[j]].x;
      vertices[j][1] = mesh->mVertices[face.mIndices[j]].y;
      vertices[j][2] = mesh->mVertices[face.mIndices[j]].z;
    }
    scene_triangle->triangle = Triangle(vertices[0], vertices[1], vertices[2]);
    // TODO remove
    Vec3 normals[3];
    for (int j = 0; j < 3; ++j) {
      normals[j][0] = mesh->mNormals[face.mIndices[j]].x;
      normals[j][1] = mesh->mNormals[face.mIndices[j]].y;
      normals[j][2] = mesh->mNormals[face.mIndices[j]].z;
      scene_triangle->normals[j] = normals[j];
    }
    scene_triangle->material = model.materials[mesh->mMaterialIndex];
    model.scene_triangles.push_back(scene_triangle);
  }
}

}  // namespace

bool loadModel(const std::string& file, Model& out_model,
               NormalType normal_type) {
  ModelLoader loader;
  if (!loader.load(file, normal_type)) {
    return false;
  }
  out_model = std::move(loader.model);
  return true;
}
