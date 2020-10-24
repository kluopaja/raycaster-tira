#include <iostream>

#include "renderer.h"

const std::string kHelp =
"Instructions:\n"
"./binary [input_file] [--verbose]";
int main(int argc, char* argv[]) {
  if(argc <= 1) {
    std::cout << kHelp << std::endl;
  }
  bool verbose = 0;
  if(argc >= 3) {
    if(std::string(argv[2]) == "--verbose") {
      verbose = 1;
    }
  }
  std::string file_name(argv[1]);
  executeInputFile(file_name, verbose);

  // this image can be used to check if your image viewer can
  // display the PPM files correctly
  if(verbose) {
    generateGammaTestImage(1000, 1000).savePPM("gamma_test.ppm");
  }
}
