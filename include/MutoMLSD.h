/* Use Maximum Likelihood of Scattering and Displacement (MLSD or ML/EM) method to get muon image
    + load initial image e.g. PoCA, or random generated noise;
    + for all muons, estimate Wij based on scattering theory and Siddon voxel path;
    + summing weight for all the recorded muons; 

   reference: 
    Schultz, L. J. O. E. (2003). COSMIC RAY MUON RADIOGRAPHY. PhD thesis
    Schultz, L. J., Blanpied, G. S., Borozdin, K. N., Fraser, A. M., Hengartner, N. W., Klimenko, A. V, … Sossong, M. J. (2007). 
    Statistical reconstruction for cosmic ray muon tomography. IEEE Transactions on Image Processing, 16(8), 1985–1993. 
    https://doi.org/10.1109/TIP.2007.901239

    Some implementation details
    + PoCA path is calculated once before MLSD iteration. Weight matrix are stored in 2-D vector.
    + for memory concern:
      (1 million muons and 100x100x100 grid will approx. consume 4 GB of memory if using double)
    + median is implemented by std::mth_element, while averaging some percentile still require sorting vectors
*/

#pragma once

#include "MutoTypes.h"
#include "MutoVReconstruction.h"
#include "MutoConfig.h"
#include "MutoSiddon.h"
#include <chrono>

class MutoMLSD : public MutoVReconstruction {

public:
    // default constructor, set default values for MLSD process
    MutoMLSD(); 
    // construct with json configuration
    MutoMLSD(json); 

    virtual ~MutoMLSD();

    // implementation of reconstruction function usin MLSD method
    virtual Image reconstruct(const MutoMuonData &);

    // re-configure the MLSD method
    void configure(json);
    VoxelGrid getCurrentGrid() { return fGrid; } // temporary for testing

private:
    // private members 
    VoxelGrid fGrid;
    
    bool fUpdateWithMean;
    int fPercentileMean;

    bool fStraightLine;
    bool fEarlyStop;
    bool fHasEnergy;
    bool fLogPrint;
    
    int fNTotal;
    int fNScatInGrid;
    int fNIteration;

    int fInitMethod;
    MTfloat fInitValue;
    std::string fInitImageFile;
    
    MTfloat fAngMag;
    MTfloat fEPS;

    std::chrono::steady_clock::time_point fTimeStart; // use steady clock

    // initialization
    Image initailizeImage();
    // convert (x, y, z) index to actual index in memory 
    inline MTindex xyz2index(MTindex x, MTindex y, MTindex z) { return x + y*fGrid.nx + z*fGrid.nx*fGrid.ny; }

    void gatherInformation(const MutoMuonData &, std::vector<MTindex>&, std::vector<MTfloat>&, std::vector<MTfloat>&, 
                           std::vector<MTfloat>&, std::vector<MTfloat>&, std::vector<MTfloat>&, std::vector<Vector3>& );
    Vector3 getPointOfClosestApproach(const RayData&);
    std::pair<MTfloat, MTfloat> getScatteringAngle(const Vector3&, const Vector3&);
    std::pair<MTfloat, MTfloat> getDisplacement(const RayData&);

    // functions for weight initialization
    //
    // get all paths and tail lengths from muons
    MTindex getMuonPathAndTail(const RayData&, const Vector3&, std::vector<VoxelData>&, std::vector<MTfloat>&, MutoSiddon&);
    // generate weight matrix from ONE muon (index i)
    MTindex generateWeightMatrix(std::vector<std::pair<MTindex, Mat2x2>>&, const std::vector<VoxelData>&, const std::vector<MTfloat>&);
    // allocate memory for S_values matrix
    void allocateSValueMatrix(std::vector<std::vector<MTfloat>>&, const std::vector<VoxelData>&);

    // functions for the update of MLSD iterations
    //
    // get covariance matrix for ONE muon (index i)
    Mat2x2 getCovarianceMatrix(const Image&, const std::vector<std::pair<MTindex, Mat2x2>>&, MTfloat);
    // update S matrix for ONE muon (for index i, update all j such that Lij != 0)
    // return the cost function value of current muon
    MTfloat updateSValueMatrix(std::vector<std::vector<MTfloat>>&, std::vector<MTindex>& ,const Image&, const Mat2x2&, 
                            const std::vector<std::pair<MTindex, Mat2x2>>&, MTfloat, MTfloat, MTfloat, MTfloat, MTfloat);
    // update one voxel of image by Sj
    MTfloat updateImageVoxel(std::vector<MTfloat>&, MTindex);

    // functions for log information
    void logMessage(const std::string&, bool=false); // second parameter check for timing 
    void logResults();
};