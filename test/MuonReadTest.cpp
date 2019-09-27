#include <iostream>
#include <chrono>
#include <random>

#include "MutoTypes.h"
#include "MutoFile.h"

int main () {
    int nPrint = 10;

    std::vector<std::string> files;
    files.push_back("/data/TUMUTY/data/old_processed/20180410_AlPbPowder_10d/raw_data/poscalieng180411.bin");
    std::vector<MTfloat> zlayers {2941.0, 2485.0, 2102.0, 1194.0, 584.0, 163.0};

    // start timing
    auto ts = std::chrono::steady_clock::now();

    MutoMuonData data = MutoFile::getMuonData(files, zlayers, 3);

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << "Reading file took: " << time_span.count() << " seconds. " << std::endl;

    for (int i=0; i<nPrint; ++i) {
        RayData ray = data[i];
        std::cout << ray.pin << " " << ray.din << std::endl;
        std::cout << ray.pout << " " << ray.dout << std::endl;

    }

    return 0;
}