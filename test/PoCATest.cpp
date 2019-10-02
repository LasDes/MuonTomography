#include <iostream>
#include <chrono>
#include <random>

#include "MutoTypes.h"
#include "MutoFile.h"
#include "MutoPoCA.h"

int main (int argc, char** argv) {
    MutoConfig config;

    if (argc == 1) {
        config.loadFile("example.json");
    } else {
        config.loadFile(argv[1]);
    }

    std::vector<std::string> files = config.get("data_files");
    std::vector<MTfloat> zlayers = config.get("detector_height");
    int nMiddle = config.get("detector_middle");

    // start timing
    auto ts = std::chrono::steady_clock::now();

    MutoMuonData data = MutoFile::getMuonData(files, zlayers, nMiddle);

    MutoPoCA poca(config.get("poca"));
    Image img = poca.reconstruct(data);

    std::string output = config.get("output_file");

    ImageHeader header;
    header.grid = poca.getCurrentGrid();
    header.fPrecision = 4;
    header.method = "poca";
    header.nIter = 0;
    header.nRay = data.size();

    MutoFile::saveImage(img, header, output);

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << "PoCA took: " << time_span.count() << " seconds. " << std::endl;

    return 0;
}