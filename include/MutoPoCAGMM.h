/* PoCA reconstruction using Gaussian Mixture Model to fit
    + record each PoCA point values in voxels 
    + at the end of reconstruction, fit scattering angles in each voxel to GMM
*/

#pragma once

#include "MutoTypes.h"
#include "MutoPoCA.h"
#include "MutoConfig.h"

class MutoPoCAGMM : public MutoPoCA {

public:
    // default constructor, set default values for PoCA process
    MutoPoCAGMM(); 
    // construct with json configuration
    MutoPoCAGMM(json); 

    virtual ~MutoPoCAGMM();

    // implementation of reconstruction function usin PoCA method
    virtual Image reconstruct(const MutoMuonData &);

private:
    // private members 
    std::vector<MTfloat> fKList; // list of multipliers for each energy group
    double fTolerance; // tolerance for iterations
    int fMaxIteration; // max iteration steps

    // Maximum likelihood fit to GMM
    MTfloat fitGMM(std::vector<MTfloat> data, MTfloat l_init);
};