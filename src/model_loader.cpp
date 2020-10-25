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
  unsigned int normal_flag = 0;
  switch (normal_type) {
    case NormalType::kSmooth:
      normal_flag = aiProcess_GenSmoothNormals;
      break;
    case NormalType::kRough:
      normal_flag = aiProcess_GenNormals;
      break;
    default:
      assert(0);
  }
  if (normal_type == NormalType::kSmooth) {
  }
  const aiScene* ai_scene = importer.ReadFile(
      file,
      aiProcess_Triangulate |  // every face will have <= 3 vertices
          normal_flag | aiProcess_ValidateDataStructure |
          aiProcess_FindDegenerates);  // finds degenerate things (e.g.
                                       // triangles)
  if (!ai_scene) {
    std::cerr << "ERROR while loading the model: " << importer.GetErrorString()
              << std::endl;
    std::exit(0);
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
  Material new_material;
  aiColor3D diffuse(0.0, 0.0, 0.0);
  aiColor3D emitted(0.0, 0.0, 0.0);
  aiColor3D specular(0.0, 0.0, 0.0);
  float specular_exp;
  // from http://assimp.sourceforge.net/lib_html/materials.html:
  // "Scales the specular color of the material.
  //  This value is kept separate from the specular color by most modelers,
  //   and so do we."
  //
  // So this is just used to scale the 'specular' color vector
  float specular_strength;
  aiColor3D transparent(0.0, 0.0, 0.0);
  // I assume this is similar to specular_strength in that the
  // 'transparent' needs to be scaled with '1-opacity'
  float opacity;
  float index_of_refraction;
  ai_material->Get(AI_MATKEY_COLOR_DIFFUSE, diffuse);
  ai_material->Get(AI_MATKEY_COLOR_AMBIENT, emitted);
  ai_material->Get(AI_MATKEY_COLOR_SPECULAR, specular);
  ai_material->Get(AI_MATKEY_SHININESS, specular_exp);
  ai_material->Get(AI_MATKEY_SHININESS_STRENGTH, specular_strength);
  ai_material->Get(AI_MATKEY_COLOR_TRANSPARENT, transparent);
  ai_material->Get(AI_MATKEY_OPACITY, opacity);
  ai_material->Get(AI_MATKEY_REFRACTI, index_of_refraction);
  new_material.diffuse = Vec3(diffuse.r, diffuse.g, diffuse.b);
  new_material.emitted = Vec3(emitted.r, emitted.g, emitted.b);
  // NOTE for some reason the specular_strength is often 0 (???)
  // therefore will be ignored for now
  // new_material.specular = specular_strength * Vec3(specular.r, specular.g,
  //                                                 specular.b);
  new_material.specular = Vec3(specular.r, specular.g, specular.b);
  new_material.specular_exp = specular_exp;
  new_material.transparent =
      (1.0f - opacity) * Vec3(transparent.r, transparent.g, transparent.b);
  new_material.index_of_refraction = index_of_refraction;
  model.materials.pushBack(new_material);
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
    SceneTriangle scene_triangle;
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
    scene_triangle.triangle = Triangle(vertices[0], vertices[1], vertices[2]);
    // TODO remove
    Vec3 normals[3];
    for (int j = 0; j < 3; ++j) {
      normals[j][0] = mesh->mNormals[face.mIndices[j]].x;
      normals[j][1] = mesh->mNormals[face.mIndices[j]].y;
      normals[j][2] = mesh->mNormals[face.mIndices[j]].z;
      scene_triangle.normals[j] = normals[j];
    }
    scene_triangle.material = model.materials[mesh->mMaterialIndex];
    model.scene_triangles.pushBack(scene_triangle);
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
