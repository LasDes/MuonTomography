/* Use Point of Closest Approach (PoCA) method to get muon image
    + get the closest point btw the in/out ray as scattering point;
    + distribute weight to each passing voxel;
    + summing weight for all the recorded muons; 

   reference: 
    Schultz, L. J. O. E. (2003). COSMIC RAY MUON RADIOGRAPHY. PhD thesis
    Online resource about distance btw lines: http://geomalgorithms.com/a07-_distance.html
*/

#pragma once

#include "MutoTypes.h"
#include "MutoVReconstruction.h"
#include "MutoConfig.h"

class MutoPoCA : public MutoVReconstruction {

public:
    // default constructor, set default values for PoCA process
    MutoPoCA(); 
    // construct with json configuration
    MutoPoCA(json); 

    virtual ~MutoPoCA();

    // implementation of reconstruction function usin PoCA method
    virtual Image reconstruct(const MutoMuonData &);

    // re-configure the PoCA method
    VoxelGrid getCurrentGrid() { return fGrid; } // temporary for testing

private:
    // private members 
    VoxelGrid fGrid;
    bool fLengthWeighted;
    bool fStraightLine;
    bool fHasEnergy;
    bool fLogPrint;
    int fNTotal;
    int fNNotParallel;
    int fNScatInGrid;
    MTfloat fAngMag;
    MTfloat fEPS;

    // some helper functions
    Image initailizeImage();
    Vector3 getPointOfClosestApproach(const RayData&);
    MTfloat getScatAngleSquare(const Vector3&, const Vector3&);
    void logResults();
};