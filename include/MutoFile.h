/* Provide file read and write method for muon data and image data */

#include <iostream>
#include <fstream>

#include "MutoTypes.h"

class MutoFile {
public:
    MutoFile() {}
    ~MutoFile() {}

    // muon ray data file handling (process data)
    static MutoMuonData getMuonData(std::string filename, std::vector<MTfloat> zlayers, MTindex nMiddle);
    // writting image data to file
    static bool saveImage(const Image& img, const ImageHeader& header, std::string filename);

    static std::string getCurrentPath();
private:
    static std::string fCwd;
};