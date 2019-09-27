#include <iostream>
#include <chrono>
#include <random>

#include "MutoTypes.h"
#include "MutoFile.h"

int main () {
    int nPrint = 100;

    // std::string file = "/data/CoTAL/CoTALMuon/MLSDEM/build/CoTALCosmicMLSDEM_simu_20mm_50n.raw"; // use past data as initial input
    std::string file = "testimage_double.bin"; 

    // start timing
    auto ts = std::chrono::steady_clock::now();

    auto image = MutoFile::loadImage(file);

    Image img = image.first;
    ImageHeader header = image.second;
    std::cout << header.fPrecision << std::endl;

    // bool fSave = MutoFile::saveImage(img, header, "testimage_double.bin");

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << "Reading and writing image took: " << time_span.count() << " seconds. " << std::endl;

    for (int i=0; i<nPrint; ++i) {
        std::cout << image.first(i) << " " << image.first(i+image.second.grid.nx) << std::endl;
    }

    return 0;
}