/* Straright Line Reconstruction (SLR) method to get muon image
    + get incoming and outgoing points of voxel region
    + collect voxels along the straight line btw points using Siddon algorithm
    + assign mean scattering angle to each voxel in the collection
*/

#pragma once

#include "MutoTypes.h"
#include "MutoVReconstruction.h"
#include "MutoConfig.h"

class MutoSLR : public MutoVReconstruction {

public:
    // default constructor, set default values for PoCA process
    MutoSLR(); 
    // construct with json configuration
    MutoSLR(json); 

    virtual ~MutoSLR();

    // implementation of reconstruction function usin SLR method
    virtual Image reconstruct(const MutoMuonData &);

    // re-configure the PoCA method
    VoxelGrid getCurrentGrid() { return fGrid; } // temporary for testing

protected:
    // private members 
    VoxelGrid fGrid;
    bool fLengthWeighted;
    bool fHasEnergy;
    bool fLogPrint;
    int fNTotal;
    MTfloat fAngMag;
    MTfloat fEPS;

    // some helper functions
    Image initailizeImage();
    void logResults();
};