#include "MutoMAP.h"
#include "MutoFile.h"
#include "MutoSiddon.h"

// --------------------------------------------------------------------
// Constructors and destructor
// --------------------------------------------------------------------
MutoMAP::MutoMAP() : MutoMLSD(){
    MTfloat fBeta = 1.0;
    MTfloat fP = 1.0;

    // some flags to enable addition features
    bool fUseAllNeighbours = false;
    bool fAddObjectPrior = false;
}

MutoMAP::MutoMAP(json config) : MutoMLSD(config) {
    fBeta = config.value("beta", 1.0);
    fP = config.value("p", 1.0);
    fAddObjectPrior = config.value("add_object_prior", false);
    fUseAllNeighbours = config.value("all_neighbours", false);
}

MutoMAP::~MutoMAP() {}

// --------------------------------------------------------------------
// Core interface of MLSD method
// --------------------------------------------------------------------
Image MutoMAP::reconstruct(const MutoMuonData& rays) {
    logMessage("==================================================");
    logMessage("|| MLSD reconstruction start                    ||");
    logMessage("==================================================");

    // start timing
    fTimeStart = std::chrono::steady_clock::now();

    Image img = initailizeImage(); // image before any ML/EM iteration
    fNScatInGrid = 0;
    fNTotal = 0 < fNTotal ? fNTotal : rays.size();

    logMessage("Initial image was loaded / generated." , true);

    // ------------------------------------------
    // Prepare scattering & displacement data
    // ------------------------------------------

    // some vectors containing basic information from rays
    std::vector<MTindex> indexInWorld; // index of rays scattering inside grid system
    std::vector<MTfloat> xScatterAngle; // scattering angle in X-Z plane
    std::vector<MTfloat> yScatterAngle; // scattering angle in Y-Z plane
    std::vector<MTfloat> xDisplacement; // displacement in X-Z plane
    std::vector<MTfloat> yDisplacement; // displacement in Y-Z plane
    std::vector<MTfloat> muonMomemtum; // momentum of muon if energy information is available
    std::vector<Vector3> scatPoint; // scattering points (may consume more memory)

    // loop once to get basic information
    gatherInformation(rays, indexInWorld, xScatterAngle, xDisplacement, yScatterAngle, yDisplacement, muonMomemtum, scatPoint);

    logMessage("Muon scattering data was generated." , true);

    // prepare memory space for intermediate results
    std::vector<VoxelData> path; // path of each muon going through the voxels
    std::vector<MTfloat> tail; // length of path from current voxel to the exit

    int nReserved = static_cast<int>(fGrid.nz * 1.5); // number of element reserved for efficiency
    path.reserve(nReserved);
    tail.reserve(nReserved);

    // weight matrix of each ray and voxel index (this is a large sparse matrix, but will save time looping through rays)
    std::vector<std::vector<std::pair<MTindex, Mat2x2>>> W_matrix; 
    W_matrix.resize(indexInWorld.size()); // resize outer dimension to the NUMBER OF RAYS
    for (auto &W : W_matrix) { W.reserve(nReserved); } // preallocate memory for inter dimension

    // Initialize matrix of S_value in memory to prevent frequent memory operation in MLSD iteration 
    std::vector<std::vector<MTfloat>> S_value; // value of Sij
    S_value.resize(img.size()); // resize outer dimension to enable preallocation of 2-D vectors
    for (auto &v : S_value) { v.reserve(nReserved); } // this may consume quite a few memory
    
    MutoSiddon siddon(fGrid);

    // loop through muon rays       
    // generate weight matrix and store it in memory 
    for (MTindex i=0; i<indexInWorld.size(); ++i) {
        // update ray path and tail length
        MTindex nHitVoxels = getMuonPathAndTail(rays[indexInWorld[i]], scatPoint[i], path, tail, siddon);
        // generate weight matrix
        MTindex nWMatrix = generateWeightMatrix(W_matrix[i], path, tail);      
        // allocate memory for S_value
        allocateSValueMatrix(S_value, path);
    }

    logMessage("Weight matrix was generated." , true);

    // ------------------------------------------
    // MLSD iteration starts here
    // ------------------------------------------

    // prepare data container for MLSD iterations
    Image imgLast = img; // image that store results from last iteration

    Mat2x2 Sigma, Sigma_inv, ErrorM; // covariance matrix 
    MTfloat Sigma_det, costFuncVal, costFuncValLast; //  determinant of cov. matrix and cost function per iteration

    // ! NOTE: THIS MATRIX IS IMPORTANT FOR REAL MEASUREMENT
    // matrix accounting for measurement errors 
    ErrorM(0,0) = 4.0*fErrorXY*fErrorXY / (fErrorZo*fErrorZo); // Ep = 1mm, dz0 = dz1 = 1000mm
    ErrorM(1,1) = 2.0*fErrorXY*fErrorXY * (1.0 + (fErrorZi/fErrorZo+1.0) * fErrorZi/fErrorZo);
    ErrorM(0,1) = ErrorM(1,0) = 2.0*fErrorXY*fErrorXY * fErrorZi/(fErrorZo*fErrorZo);

    // matrix store current index of each voxel in S_value matrix
    std::vector<MTindex> S_index;
    S_index.resize(img.size());

    logMessage("\nStart EM iteration to find the Maximum Likelihood ... ");

    // start iteration
    for (MTindex iter=0; iter<fNIteration; ++iter) {
        // clear index of S matrix
        std::fill(S_index.begin(), S_index.end(), 0);

        // reset cost function value
        costFuncValLast = costFuncVal;
        costFuncVal = 0.0;

        // update covariance matrix Sigma by loop through Wi
        for (MTindex i=0; i<indexInWorld.size(); ++i) {
            MTfloat pr = fHasEnergy ? muonMomemtum[i] : muonMomemtum[0];
            Sigma = getCovarianceMatrix(img, W_matrix[i], pr);
            Sigma += ErrorM;
            // update Sij for this muon return the value of cost function
            costFuncVal += updateSValueMatrix(S_value, S_index, img, Sigma, W_matrix[i], pr, 
                                              xScatterAngle[i], yScatterAngle[i], xDisplacement[i], yDisplacement[i]);
        }
        // update image using mean (with percentile) or median value of SMatrix
    #pragma omp parallel for
        for (MTindex j=0; j<img.size(); ++j) {
            img(j) = updateImageVoxel(S_value[j], S_index[j]);
        }
        // record cost function values
        costFuncVal /= static_cast<MTfloat>(indexInWorld.size());

        // log for current iteration
        logMessage("Iteration #" + std::to_string(iter) + " finished. ", true);
        logMessage("\tcost function = " + std::to_string(costFuncVal));

        // early stop if cost function is not decreasing 
        if (fEarlyStop && costFuncValLast < costFuncVal) {
            logMessage("Early stop at iteration#" + std::to_string(iter) + ". Minimum cost function = " + std::to_string(costFuncValLast));
            break;
        }
    }

    if (fLogPrint) {
        logResults();
    }

    return img;
}