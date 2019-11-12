#include <iostream>
#include <chrono>
#include <random>

#include "MutoTypes.h"
#include "MutoFile.h"
#include "MutoPoCA.h"
#include "MutoPoCut.h"
#include "MutoMLSD.h"
#include "MutoMAP.h"
#include "MutoSLR.h"
#include "MutoSLRIter.h"


int main (int argc, char** argv) {
    // load JSON configuration file from user input
    MutoConfig config;
    if (argc == 1) {
        config.loadFile("example.json");
    } else {
        config.loadFile(argv[1]);
    }

    // start timing
    auto ts = std::chrono::steady_clock::now();

    // load data into memory
    std::vector<std::string> files = config.get("data_files");
    std::vector<MTfloat> zlayers = config.get("detector_height");
    int nMiddle = config.get("detector_middle");
    MutoMuonData data = MutoFile::getMuonData(files, zlayers, nMiddle);

    MutoVReconstruction* rec;
    std::string method;
    ImageHeader header;

    if (config.exist("poca")) { 
        method = "poca"; 
        rec = new MutoPoCA(config.get("poca"));   
        header.grid = ((MutoPoCA*) rec) -> getCurrentGrid();
        header.nIter = 1;
    } else if (config.exist("mlsd")) {
        method = "mlsd"; 
        rec = new MutoMLSD(config.get("mlsd"));
        header.grid = ((MutoMLSD*) rec) -> getCurrentGrid();
        header.nIter = config.get("mlsd")["num_iteration"];
    } else if (config.exist("map")) {
        method = "map"; 
        rec = new MutoMAP(config.get("map"));
        header.grid = ((MutoMAP*) rec) -> getCurrentGrid();
        header.nIter = config.get("map")["num_iteration"];
    } else if (config.exist("slr")) { 
        method = "slr"; 
        rec = new MutoSLR(config.get("slr"));   
        header.grid = ((MutoSLR*) rec) -> getCurrentGrid();
        header.nIter = 1;
    } else if (config.exist("slr_iter")) { 
        method = "slr_iter"; 
        rec = new MutoSLRIter(config.get("slr_iter"));   
        header.grid = ((MutoSLRIter*) rec) -> getCurrentGrid();
        header.nIter = 1;
    } else if (config.exist("pocut")) { 
        method = "pocut"; 
        rec = new MutoPoCut(config.get("pocut"));   
        header.grid = ((MutoPoCut*) rec) -> getCurrentGrid();
        header.nIter = 1;
    } 

    Image img = rec->reconstruct(data);

    std::string output = config.get("output_file");
    
    header.fPrecision = 4;
    header.method = method;
    header.nRay = data.size();

    MutoFile::saveImage(img, header, output);

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << "Reconstruction with <" << method << "> took: " << time_span.count() << " seconds. " << std::endl;

    return 0;
}