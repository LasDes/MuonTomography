/* PoCA reconstruction with threshold of scattering angles
    + record each PoCA point values in voxels 
    + set percentile of maximum allowed scattering angle (weighted)
    + calculate statistics of scattering angles below threshold for each voxel 
*/

#pragma once

#include "MutoTypes.h"
#include "MutoPoCA.h"
#include "MutoConfig.h"

class MutoPoCut : public MutoPoCA {

public:
    // default constructor, set default values for PoCA process
    MutoPoCut(); 
    // construct with json configuration
    MutoPoCut(json); 

    virtual ~MutoPoCut();

    // implementation of reconstruction function usin PoCA method
    virtual Image reconstruct(const MutoMuonData &);

private:
    // private members 
    double fPercentile;
};