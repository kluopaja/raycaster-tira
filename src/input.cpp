#include "input.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <optional>
namespace {

const std::string kCommentString = "//";
// Enums corresponding different input file parts
enum class Section {
  kCamera,
  kModel,
  kPointLight,
  kEnvironmentLight,
  kRender,
  kImage,
  N_ENUM_VALUES,
};
enum class CameraField {
  kPos,
  kFront,
  kUp,
  kXFov,
  kYFov,
  N_ENUM_VALUES,
};
enum class ModelField {
  kFile,
  kPos,
  kNormal,
  N_ENUM_VALUES,
};
enum class PointLightField {
  kPos,
  kColor,
  N_ENUM_VALUES,
};
enum class EnvironmentLightField {
  kColor,
  kType,
  kDirection,
  kExp,
  N_ENUM_VALUES,
};
enum class RenderField {
  kWidth,
  kHeight,
  kPixelRays,
  kDepth,
  kBranching,
  N_ENUM_VALUES,
};
enum class ImageField {
  kFile,
  kTruncate,
  kScaleMax,
  N_ENUM_VALUES,
};

// Functions for applying the actual mapping between enums and strings
// TODO also overload operator>> for this?
std::string enumToStr(Section s) {
  std::string section_names[(std::size_t)Section::N_ENUM_VALUES];
  section_names[(std::size_t)Section::kCamera] = "camera";
  section_names[(std::size_t)Section::kModel] = "model";
  section_names[(std::size_t)Section::kPointLight] = "point_light";
  section_names[(std::size_t)Section::kEnvironmentLight] = "environment_light";
  section_names[(std::size_t)Section::kRender] = "render";
  section_names[(std::size_t)Section::kImage] = "image";
  return section_names[(std::size_t)s];
}
std::optional<Section> strToSection(const std::string& s) {
  for(std::size_t i = 0; i < (std::size_t)Section::N_ENUM_VALUES; ++i) {
    if (s == enumToStr(static_cast<Section>(i))) {
      return static_cast<Section>(i);
    }
  }
  return std::nullopt;
}
std::string enumToStr(CameraField x) {
  std::string camera_names[(std::size_t)CameraField::N_ENUM_VALUES];
  camera_names[(std::size_t)CameraField::kPos] = "pos";
  camera_names[(std::size_t)CameraField::kFront] = "front";
  camera_names[(std::size_t)CameraField::kUp] = "up";
  camera_names[(std::size_t)CameraField::kXFov] = "x_fov";
  camera_names[(std::size_t)CameraField::kYFov] = "y_fov";
  return camera_names[(std::size_t)x];
}
std::string enumToStr(ModelField x) {
  std::string model_names[(std::size_t)ModelField::N_ENUM_VALUES];
  model_names[(std::size_t)ModelField::kFile] = "file";
  model_names[(std::size_t)ModelField::kPos] = "pos";
  model_names[(std::size_t)ModelField::kNormal] = "normal";
  return model_names[(std::size_t)x];
}
std::string enumToStr(PointLightField x) {
  std::string point_light_names[(std::size_t)PointLightField::N_ENUM_VALUES];
  point_light_names[(std::size_t)PointLightField::kPos] = "pos";
  point_light_names[(std::size_t)PointLightField::kColor] = "color";
  return point_light_names[(std::size_t)x];
}
std::string enumToStr(EnvironmentLightField x) {
  std::string environment_light_names[(std::size_t)EnvironmentLightField::N_ENUM_VALUES];
  environment_light_names[(std::size_t)EnvironmentLightField::kColor] = "color";
  environment_light_names[(std::size_t)EnvironmentLightField::kType] = "type";
  environment_light_names[(std::size_t)EnvironmentLightField::kDirection] = "direction";
  environment_light_names[(std::size_t)EnvironmentLightField::kExp] = "exp";
  return environment_light_names[(std::size_t)x];
}
std::string enumToStr(RenderField x) {
  std::string render_names[(std::size_t)RenderField::N_ENUM_VALUES];
  render_names[(std::size_t)RenderField::kWidth] = "width";
  render_names[(std::size_t)RenderField::kHeight] = "height";
  render_names[(std::size_t)RenderField::kPixelRays] = "pixel_rays";
  render_names[(std::size_t)RenderField::kDepth] = "depth";
  render_names[(std::size_t)RenderField::kBranching] = "branching";
  return render_names[(std::size_t)x];
}
std::string enumToStr(ImageField x) {
  std::string image_names[(std::size_t)ImageField::N_ENUM_VALUES];
  image_names[(std::size_t)ImageField::kFile] = "file";
  image_names[(std::size_t)ImageField::kTruncate] = "truncate";
  image_names[(std::size_t)ImageField::kScaleMax] = "scale_max";
  return image_names[(std::size_t)x];
}
// Defines the mapping between NormalType values and input file strings
// Note that this is here and not where NormalType is defined
// because this mapping is limited to the file format.
std::string enumToStr(NormalType x) {
  std::string normal_type_names[(std::size_t)NormalType::N_ENUM_VALUES];
  normal_type_names[(std::size_t)NormalType::kSmooth] = "smooth";
  normal_type_names[(std::size_t)NormalType::kRough] = "rough";
  return normal_type_names[(std::size_t)x];
}
// Reading the field values is done with the overloaded operator>>
std::istream& operator>>(std::istream& in, NormalType& out) {
  std::string tmp;
  in >> tmp;
  bool found_match = 0;
  for(std::size_t i = 0; i < (std::size_t)NormalType::N_ENUM_VALUES; ++i) {
    if(tmp == enumToStr(static_cast<NormalType>(i))) {
      out = static_cast<NormalType>(i);
      found_match = 1;
    }
  }
  if(!found_match) {
    in.setstate(std::ios::failbit);
  }
  return in;
}
std::ostream& operator<<(std::ostream& out, const NormalType& a) {
  out << enumToStr(a);
  return out;
}
// Defines the mapping between EnvironmentLightType values and input
// file strings
std::string enumToStr(EnvironmentLightType x) {
  std::string environment_light_type_names[(std::size_t)EnvironmentLightType::N_ENUM_VALUES];
  environment_light_type_names[(std::size_t)EnvironmentLightType::kUniform] = "uniform";
  environment_light_type_names[(std::size_t)EnvironmentLightType::kDirected] = "directed";
  return environment_light_type_names[(std::size_t)x];
}
std::istream& operator>>(std::istream& in, EnvironmentLightType& out) {
  std::string tmp;
  in >> tmp;
  bool found_match = 0;
  for(std::size_t i = 0; i < (std::size_t)EnvironmentLightType::N_ENUM_VALUES; ++i) {
    if(tmp == enumToStr(static_cast<EnvironmentLightType>(i))) {
      out = static_cast<EnvironmentLightType>(i);
      found_match = 1;
    }
  }
  if(!found_match) {
    in.setstate(std::ios::failbit);
  }
  return in;
}
std::ostream& operator<<(std::ostream& out, const EnvironmentLightType& a) {
  out << enumToStr(a);
  return out;
}

// Reports an error for too many definitions of some section
void printMultiSectionError(Section section, int line) {
  std::cerr << "Line " << line
            << ": Second definition of ";
  std::cerr << enumToStr(section);
  std::cerr << ". Only one is allowed per file" << std::endl;
}
// Reports an error for many identically named fields within a section
void printMultiFieldError(const std::string& name, int line) {
  std::cerr << "Line " << line << ": Second definition of field ";
  std::cerr << name << ". Only one is allowed per section." << std::endl;
}
// Reports an error for incorrectly formatted field value
void printFieldParseError(const std::string& name, int line) {
  std::cerr << "Line " << line << ": Failed parsing value of field ";
  std::cerr << name << "." << std::endl;
}
// Stores one line of the input. `parsed` is set to 1 after
// the line has been processed
struct InputLine {
  bool parsed = 0;
  std::string content;
};
// Reads the `file_name` file to a Vector of strings each 
// element of the Vector corresponding to one line in the file
Vector<InputLine> readLines(const std::string& file_name) {
  std::ifstream fin(file_name);
  assert(fin.good() && "Opening the input file failed");
  Vector<InputLine> out;
  while (fin.good()) {
    std::string s;
    std::getline(fin, s);
    if (fin.good()) {
      InputLine tmp;
      tmp.content = s;
      out.pushBack(tmp);
    }
  }
  return out;
}
// Returns the first substring separated by whitespaces on `line`
std::string extractFirstString(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string first_part;
  ss >> first_part;
  return first_part;
}
// Uses overloaded operator>> to extract the value part of a field
// Either returns the correctly parsed type or std::nullopt
template<class T> 
std::optional<T> extractFieldValue(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string first_part;
  ss >> first_part;
  T out;
  if(ss >> out) {
    return out;
  }
  return std::nullopt;
}
// Stores a section of type `type` starting at position `position`
struct SectionBegin {
  Section type;
  std::size_t position;
};
// Determines if `line` marks a start of a section
// Returns the corresponding SectionType or std::nullopt
std::optional<Section> parseSectionBegin(const std::string& line) {
  std::string name = extractFirstString(line);
  return strToSection(name);
}
// Finds the types and starting positions of all sections in `input`
// Marks the corresponding lines as parsed
Vector<SectionBegin> parseSectionBegins(Vector<InputLine>& input) {
  Vector<SectionBegin> out;
  for(std::size_t i = 0; i < input.size(); ++i) {
    std::optional<Section> section = parseSectionBegin(input[i].content);
    if (!section) continue;
    input[i].parsed = 1;
    out.pushBack({*section, i});
  }
  return out;
}
void printSectionBegins(const Vector<SectionBegin>& section_begins) {
  std::cout << "Found the following sections from the input file:" << std::endl;
  for(const SectionBegin& x: section_begins) {
    std::cout << "\t";
    std::cout << "Section of type ";
    std::cout << std::left << std::setw(20);
    std::cout <<  enumToStr(x.type) << " ";
    std::cout << "on line: " << x.position + 1 << std::endl;
  }
}
// Checks that there are no multiple definitions
// of sections `kCamera`, `kEnvironmentLight` or `kRender`
// Exits the program if multiple definitions are found
void validateSections(const Vector<SectionBegin>& section_begins) {
  int n_cameras = 0;
  int n_environment_lights = 0;
  int n_renders = 0;
  for(auto section: section_begins) {
    switch (section.type) {
      case Section::kCamera:
        ++n_cameras;
        break;
      case Section::kEnvironmentLight:
        ++n_environment_lights;
        break;
      case Section::kRender:
        ++n_renders;
        break;
      default:
        break;
    }
    if(n_cameras > 1 || n_environment_lights > 1 || n_renders > 1) {
      printMultiSectionError(section.type, section.position);
      std::exit(0);
    }
  }
}
// Stores the location of a section
struct SectionRange {
  std::size_t begin = 0;
  std::size_t end = 0;
};
// Finds a field with name `name` from the SectionRange `range`
// Checks that only one field is found. Otherwise
// prints an error and exits the program.
std::optional<std::size_t> findFieldByName(SectionRange range,
                                           const Vector<InputLine>& input,
                                           const std::string& name) {
  Vector<std::size_t> positions;
  for(std::size_t i = range.begin; i < range.end; ++i) {
    if(extractFirstString(input[i].content) == name) {
      positions.pushBack(i);
    }
  }
  if(positions.size() == 0) {
    return std::nullopt;
  }
  if(positions.size() == 1) {
    return positions[0];
  }
  printMultiFieldError(input[positions[1]].content, positions[1]);
  std::exit(0);
}
// Extracts the value of field `name` in the `range` to `out`
// and marks the corresponding line parsed.
//
// If no matching field is found, then `out` stays intact.
//
// An out parameter is used so that its type can be deduced.
template<class T>
void extractField(SectionRange range, Vector<InputLine>& input,
                       const std::string& name, bool verbose, T& out) {
  std::optional<std::size_t> pos = findFieldByName(range, input, name);
  if(pos) {
    if (verbose) {
      std::cout << "\t\tFound field " << name << std::endl;
    }
    std::optional<T> field_value = extractFieldValue<T>(input[*pos].content);
    if(field_value) {
      input[*pos].parsed = 1;
      out = *field_value;
    }
    else {
      printFieldParseError(name, *pos);
      std::exit(0);
    }
  }
  else {
    if(verbose) {
      std::cout << "\t\tDidn't find field " << name 
                << ". Using default value." << std::endl;
    }
  }
}
// Functions for parsing the different sections.
CameraConfig parseCamera(SectionRange range, Vector<InputLine>& input, bool verbose) {
  CameraConfig out;
  extractField(range, input, enumToStr(CameraField::kPos), verbose, out.pos);
  extractField(range, input, enumToStr(CameraField::kFront), verbose, out.front);
  extractField(range, input, enumToStr(CameraField::kUp), verbose, out.up);
  // TODO replace with something!
  // Now converts the read degrees to radians for CameraConfig
  double x_fov_degrees = out.x_fov / kPi * 180.0;
  double y_fov_degrees = out.y_fov / kPi * 180.0;
  extractField(range, input, enumToStr(CameraField::kXFov), verbose, x_fov_degrees);
  extractField(range, input, enumToStr(CameraField::kYFov), verbose, y_fov_degrees);
  out.x_fov = x_fov_degrees * kPi / 180.0;
  out.y_fov = y_fov_degrees * kPi / 180.0;
  return out;
}
void validateModel(SectionRange range, const ModelConfig& config) {
  if(config.file.empty()) {
    std::cerr << "In section from line " << range.begin << " to line "
              << range.end << ": no filename found or the filename is invalid. "
              << std::endl;
    std::exit(0);
  }
}
ModelConfig parseModel(SectionRange range, Vector<InputLine>& input, bool verbose) {
  ModelConfig out;
  extractField(range, input, enumToStr(ModelField::kFile), verbose, out.file);
  extractField(range, input, enumToStr(ModelField::kPos), verbose, out.pos);
  extractField(range, input, enumToStr(ModelField::kNormal), verbose, out.normal);
  validateModel(range, out);
  return out;
}
PointLightConfig parsePointLight(SectionRange range, Vector<InputLine>& input, bool verbose) {
  PointLightConfig out;
  extractField(range, input, enumToStr(PointLightField::kPos), verbose, out.pos);
  extractField(range, input, enumToStr(PointLightField::kColor), verbose, out.color);
  return out;
}
EnvironmentLightConfig parseEnvironmentLight(SectionRange range,
                                             Vector<InputLine>& input, bool verbose) {
  EnvironmentLightConfig out;
  extractField(range, input, enumToStr(EnvironmentLightField::kColor),
               verbose, out.color);
  extractField(range, input, enumToStr(EnvironmentLightField::kType),
               verbose, out.type);
  extractField(range, input, enumToStr(EnvironmentLightField::kDirection),
               verbose, out.direction);
  extractField(range, input, enumToStr(EnvironmentLightField::kExp),
               verbose, out.exp);
  return out;
}
RenderConfig parseRender(SectionRange range, Vector<InputLine>& input, bool verbose) {
  RenderConfig out;
  extractField(range, input, enumToStr(RenderField::kWidth), verbose, out.width);
  extractField(range, input, enumToStr(RenderField::kHeight), verbose, out.height);
  extractField(range, input, enumToStr(RenderField::kPixelRays), verbose, out.pixel_rays);
  extractField(range, input, enumToStr(RenderField::kDepth), verbose, out.depth);
  extractField(range, input, enumToStr(RenderField::kBranching), verbose, out.branching);
  return out;
}
void validateImage(SectionRange range, const ImageConfig& config) {
  if(config.file.empty()) {
    std::cerr << "In section from line " << range.begin << " to line "
              << range.end << ": no filename found or the filename is invalid. "
              << std::endl;
    std::exit(0);
  }
}
ImageConfig parseImage(SectionRange range, Vector<InputLine>& input, bool verbose) {
  ImageConfig out;
  extractField(range, input, enumToStr(ImageField::kFile), verbose, out.file);
  extractField(range, input, enumToStr(ImageField::kTruncate), verbose, out.truncate);
  extractField(range, input, enumToStr(ImageField::kScaleMax), verbose, out.scale_max);
  validateImage(range, out);
  return out;
}
// Parses a section of type `section` spanning range `range` and adds the
// result to `recipe`
void parseSection(SectionRange range, Vector<InputLine>& input,
                  Section section, bool verbose, Recipe& recipe) {
  if (verbose) {
    std::cout << "\tParsing section: " << enumToStr(section)
              << " from line " << range.begin+1 << " to line "
              << range.end << std::endl;
  }
  switch (section) {
    case Section::kCamera:
      recipe.camera = parseCamera(range, input, verbose);
      break;
    case Section::kModel:
      recipe.models.pushBack(parseModel(range, input, verbose));
      break;
    case Section::kPointLight:
      recipe.point_lights.pushBack(parsePointLight(range, input, verbose));
      break;
    case Section::kEnvironmentLight:
      recipe.environment_light = parseEnvironmentLight(range, input, verbose);
      break;
    case Section::kRender:
      recipe.render = parseRender(range, input, verbose);
      break;
    case Section::kImage:
      recipe.images.pushBack(parseImage(range, input, verbose));
      break;
    default:
      assert(0);
  }
}
// Parses sections in 'section_begins' to a Recipe
Recipe parseSections(Vector<SectionBegin> section_begins,
                     Vector<InputLine>& input, bool verbose) {
  if(verbose) {
    std::cout << "Starting to parse section" << std::endl;
  }
  Recipe recipe;
  for (std::size_t i = 0; i < section_begins.size(); ++i) {
    SectionRange section;
    section.begin = section_begins[i].position;
    section.end = 0;
    if (i+1 == section_begins.size()) {
      section.end = input.size();
    }
    else {
      section.end = section_begins[i+1].position;
    }
    parseSection(section, input, section_begins[i].type, verbose, recipe);
  }
  return recipe;
}
// checks if the line starts with the `kCommentString`
bool lineIsComment(const std::string& line) {
  std::string first = extractFirstString(line);
  return first.substr(0, 2) == kCommentString;
}
// Removes the smallest suffix starting with '/'
// if no such suffix is found, then returns an empty string.
std::string getParentPath(const std::string& path) {
  if(path.size() == 0) {
    return "";
  }
  std::size_t slash_index = path.rfind("/");
  if(slash_index == std::string::npos) {
    return "";
  }
  return path.substr(0, slash_index);
}
// Adds `input_file_path` to the parsed file paths in `recipe`
void fixFilePaths(const std::string& input_file_path, Recipe& recipe) {
  std::string parent_path = getParentPath(input_file_path);
  for(auto& model: recipe.models) {
    model.file = parent_path + "/" + model.file;
  }
  for(auto& image: recipe.images) {
    image.file = parent_path + "/" + image.file;
  }
}
// Finds comments in `input` and marks them parsed
void parseComments(Vector<InputLine>& input, bool verbose) {
  if(verbose) {
    std::cout << "Starting to parse comments" << std::endl;
  }
  for(std::size_t i = 0; i < input.size(); ++i) {
    if(!input[i].parsed && lineIsComment(input[i].content)) {
      input[i].parsed = 1;
      if(verbose) {
        std::cout << "Found comment on line " << i+1 << std::endl;
      }
    }
  }
}
// Returns true if `line` contains non-whitespace characters
// false otherwise
bool lineIsEmpty(const std::string& line) {
  std::string first = extractFirstString(line);
  return first.size() == 0;
}
// Finds empty lines (only whitespaces) in `input` and marks them parsed
void parseEmptyLines(Vector<InputLine>& input, bool verbose) {
  if(verbose) {
    std::cout << "Starting to parse empty lines" << std::endl;
  }
  for(std::size_t i = 0; i < input.size(); ++i) {
    if(!input[i].parsed && lineIsEmpty(input[i].content)) {
      input[i].parsed = 1;
      if(verbose) {
        std::cout << "Found empty line " << i+1 << std::endl;
      }
    }
  }
}
// Checks that every line in `input` has been parsed
// Otherwise prints an error and exits the program
void checkFullyParsed(Vector<InputLine>& input, bool verbose) {
  if(verbose) {
    std::cout << "Starting to check that input was fully parsed" << std::endl;
  }
  for(std::size_t i = 0; i < input.size(); ++i) {
    if(!input[i].parsed) {
      std::cout << "Failed to parse line " << i+1 << std::endl;
      std::exit(0);
    }
  }
}
// TODO implement
void validateRecipe(Recipe& recipe, bool verbose) {
  if(verbose) {
    std::cout << "Starting to validate the recipe" << std::endl;
    std::cout << "NOT FULLY IMPLEMENTED" << std::endl;
  }
}
void printRecipe(Recipe& recipe) {
  std::cout << "The following recipe was parsed from the files" << std::endl;
  std::cout << recipe << std::endl;
}
}  // namespace

std::ostream& operator<<(std::ostream& out, const CameraConfig& a) {
  out << "\t";
  out << "CameraConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "pos: "   << a.pos << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "front: " << a.front << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "up: "    << a.up << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "x_fov: " << a.x_fov << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "y_fov: " << a.y_fov << std::endl;
  out << "\t";
  out << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const ModelConfig& a) {
  out << "\t";
  out << "ModelConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "file: " << a.file << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "pos: " << a.pos << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "normal: " << a.normal << std::endl;
  out << "\t";
  out << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const PointLightConfig& a) {
  out << "\t";
  out << "PointLightConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "pos: " << a.pos << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "color: " << a.color << std::endl;
  out << "\t";
  out << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const EnvironmentLightConfig& a) {
  out << "\t";
  out << "EnvironmentLightConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "color: " << a.color << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "type: " << a.type << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "direction: " << a.direction << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "exp " << a.exp << std::endl;
  out << "\t";
  out << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const RenderConfig& a) {
  out << "\t";
  out << "RenderConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "width: " << a.width << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "height: " << a.height << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "pixel_rays: " << a.pixel_rays << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "depth: " << a.depth << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "branching: " << a.branching << std::endl;
  out << "\t";
  out << ")";
  return out;
}
std::ostream& operator<<(std::ostream& out, const ImageConfig& a) {
  out << "\t";
  out << "ImageConfig(" << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "file: " << a.file << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "truncate: " << a.truncate << std::endl;
  out << "\t\t";
  out << std::left << std::setw(12);
  out << "scale_max: " << a.scale_max << std::endl;
  out << "\t";
  out << ")";
  return out;
}

std::ostream& operator<<(std::ostream& out, const Recipe& a) {
  out << "Recipe(" << std::endl;
  out << a.camera << std::endl;
  for(auto& x: a.models) {
    out << x << std::endl;
  }
  for(auto& x: a.point_lights) {
    out << x << std::endl;
  }
  out << a.environment_light << std::endl;
  out << a.render << std::endl;
  for(auto& x: a.images) {
    out << x << std::endl;
  }
  out << ")";
  return out;
}
// Loads a rendering recipe from file `file`
//
// First reads file into `input` Vector
// After this the input is parsed and every time a
// meaning of a `input` line is stored to other variables,
// the corresponding line is marked as parsed.
//
// After this, the function checks that `input`
// only contains parsed lines, empty lines or commented lines
//
// Finally, the validity of the parsed recipe is checked 
Recipe loadRecipe(const std::string& file_name, bool verbose) {
  Vector<InputLine> input = readLines(file_name);
  // Type of the section and the position where it was found
  Vector<SectionBegin> section_begins = parseSectionBegins(input);
  if(verbose) {
    printSectionBegins(section_begins);
  }
  validateSections(section_begins);
  Recipe recipe = parseSections(section_begins, input, verbose);
  fixFilePaths(file_name, recipe);
  parseComments(input, verbose);
  parseEmptyLines(input, verbose);
  checkFullyParsed(input, verbose);
  validateRecipe(recipe, verbose);
  if(verbose) {
    printRecipe(recipe);
  }
  return recipe;
}
