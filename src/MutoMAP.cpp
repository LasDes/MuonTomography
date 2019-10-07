#include "MutoMAP.h"
#include "MutoFile.h"
#include "MutoSiddon.h"

// --------------------------------------------------------------------
// Constructors and destructor
// --------------------------------------------------------------------
MutoMAP::MutoMAP() : MutoMLSD(){
    defaultInitialize();
}

MutoMAP::MutoMAP(json config) : MutoMLSD(config) {
    defaultInitialize();
    fBeta = config.value("beta", 1.0);
    fP = config.value("p", 1.0);
    fDelta = config.value("delta", 0.1);
    fAddObjectPrior = config.value("add_object_prior", false);
    fUseAllNeighbours = config.value("all_neighbours", false);
}

MutoMAP::~MutoMAP() {}

void MutoMAP::defaultInitialize() {
    fBeta = 1.0;
    fP = 1.0;
    fDelta = 0.1;

    // some flags to enable addition features
    fUseAllNeighbours = false;
    fAddObjectPrior = false;

    // others
    fUnitConv = 1.0e7;
}

// --------------------------------------------------------------------
// Core interface of MAP (ML with regularization) method
// --------------------------------------------------------------------
Image MutoMAP::reconstruct(const MutoMuonData& rays) {
    logMessage("==================================================");
    logMessage("|| MLSD reconstruction start                    ||");
    logMessage("==================================================");

    // start timing
    fTimeStart = std::chrono::steady_clock::now();

    Image img = initailizeImage(); // image before any ML/EM iteration
    MTindex imgsize = fGrid.nx * fGrid.ny * fGrid.nz;
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
    S_value.resize(imgsize); // resize outer dimension to enable preallocation of 2-D vectors
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
    // Iteration starts here
    // ------------------------------------------

    // prepare data container for MLSD iterations
    Image imgLast = img; // image that store results from last iteration
    Image imgML = img; // intermediate image that results from ML/EM algorithm

    Mat2x2 Sigma, Sigma_inv, ErrorM; // covariance matrix 
    MTfloat Sigma_det, costFuncVal, costFuncValLast; //  determinant of cov. matrix and cost function per iteration

    // ! NOTE: THIS MATRIX IS IMPORTANT FOR REAL MEASUREMENT
    // matrix accounting for measurement errors 
    ErrorM(0,0) = 4.0*fErrorXY*fErrorXY / (fErrorZo*fErrorZo); // Ep = 1mm, dz0 = dz1 = 1000mm
    ErrorM(1,1) = 2.0*fErrorXY*fErrorXY * (1.0 + (fErrorZi/fErrorZo+1.0) * fErrorZi/fErrorZo);
    ErrorM(0,1) = ErrorM(1,0) = 2.0*fErrorXY*fErrorXY * fErrorZi/(fErrorZo*fErrorZo);

    // matrix store current index of each voxel in S_value matrix
    std::vector<MTindex> S_index;
    S_index.resize(imgsize);

    logMessage("\nStart EM iteration to find the Maximum Likelihood ... ");

    // start iteration
    for (MTindex iter=0; iter<fNIteration; ++iter) {
        // ------------------------------------------
        // MLSD/EM part of each iteration 
        // ------------------------------------------

        // clear index of S matrix
        std::fill(S_index.begin(), S_index.end(), 0);

        // reserve image from last iteration (MAP optimization is inplace of `img`)
        imgLast = img;

        // reset cost function value
        costFuncValLast = costFuncVal;
        costFuncVal = 0.0;

        // update covariance matrix Sigma by loop through Wi
        for (MTindex i=0; i<indexInWorld.size(); ++i) {
            MTfloat pr = fHasEnergy ? muonMomemtum[i] : muonMomemtum[0];
            Sigma = getCovarianceMatrix(imgLast, W_matrix[i], pr);
            Sigma += ErrorM;
            // update Sij for this muon return the value of cost function
            costFuncVal += updateSValueMatrix(S_value, S_index, imgLast, Sigma, W_matrix[i], pr, 
                                              xScatterAngle[i], yScatterAngle[i], xDisplacement[i], yDisplacement[i]);
        }
        // update image using mean (with percentile) or median value of SMatrix
    #pragma omp parallel for
        for (MTindex j=0; j<imgsize; ++j) {
            imgML(j) = updateImageVoxel(S_value[j], S_index[j]);
        }

        // ------------------------------------------
        // MAP regularization starts here
        // ------------------------------------------

        // optimize each voxel one-by-one
    #pragma omp parallel for
        for (MTindex j=0; j<imgsize; ++j) { 
            if (S_index[j] == 0) {
                img(j) = 0.0;
                continue;
            }
            // prepare coefficients for the cubic equation which represent the optimize step
            MTfloat p0, p1, p2, p3;
            p0 = -imgML(j) * fUnitConv;
            p1 = 1.0;
            p2 = p3 = 0.0; // for later accumulation

            // calculation of p2 and p3
            auto weightAndLambda = calculateCoefFromNeighbours(j, imgLast, S_index[j], p2, p3);
            // solve cubic equation: p0 + p1 x + p2 x^2 + p3 x^3 = 0
            MTfloat roots[3];
            MTindex nRoots = cubicEqnSolver(roots, p3, p2, p1, p0, fEPS);

            // choose the best solution and update image
            img(j) = chooseBestEstimation(S_index[j], imgLast(j), imgML(j), roots, nRoots, weightAndLambda);
        }

        // record cost function values
        costFuncVal /= static_cast<MTfloat>(indexInWorld.size());

        // log for current iteration
        logMessage("Iteration #" + std::to_string(iter) + " finished. ", true);
        logMessage("\tMLSD cost function = " + std::to_string(costFuncVal));

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

std::pair<std::vector<MTfloat>, std::vector<MTfloat>> MutoMAP::calculateCoefFromNeighbours(MTindex j, const Image& img, MTindex Mj, MTfloat& p2, MTfloat& p3){
    // generate weights and neighbours
    std::vector<MTfloat> wm, lm;
    if (!fUseAllNeighbours) {
        wm = std::vector<MTfloat>(6, 1.0);
        lm = getSixNeighbours(j, img);        
    } else {
        MTfloat w2 = 1.0 / std::sqrt(2.0);
        MTfloat w3 = 1.0 / std::sqrt(3.0);
        for (int i=0; i<6; ++i) { wm.push_back(1.0); }
        for (int i=0; i<12; ++i) { wm.push_back(w2); }
        for (int i=0; i<8; ++i) { wm.push_back(w3); }
        lm = getAllNeighbours(j, img);
    }
    // add each voxel value to p2 and p3
    MTfloat l_j = img(j) * fUnitConv; // conversion of unit here to avoid modifying img
    for (int i=0; i<wm.size(); ++i) {
        MTfloat t1 = wm[i] * fP * std::pow(std::pow(l_j - lm[i], 2) + fDelta*fDelta, (fP-2.0)/2.0);
        p2 += t1 * 2.0 * (l_j + lm[i]);
        p3 += t1 * 4.0;
    }
    // add object prior to p3 if configured
    if (fAddObjectPrior) {
        p3 += fP * std::pow(std::pow(l_j, 2) + fDelta*fDelta, (fP-2.0)/2.0);
    }
    
    p2 *= -fBeta / static_cast<MTfloat>(Mj);
    p3 *= fBeta / static_cast<MTfloat>(Mj);

    return std::pair<std::vector<MTfloat>, std::vector<MTfloat>>{wm, lm};
}

MTfloat MutoMAP::chooseBestEstimation(MTindex Mj, MTfloat ljn, MTfloat lj_ML, MTfloat* roots, MTindex nRoots, 
                                      const std::pair<std::vector<MTfloat>, std::vector<MTfloat>>& wl){
    // scale image values
    MTfloat ljnScaled = ljn * fUnitConv;
    MTfloat ljMLScaled = lj_ML * fUnitConv;

    std::vector<MTfloat> funcValue;
    for (MTindex i=0; i<nRoots; ++i) {
        funcValue.push_back(surrogateFunctionValue(Mj, ljnScaled, ljMLScaled, roots[i], wl.first, wl.second));
    }
    funcValue.push_back(surrogateFunctionValue(Mj, ljnScaled, ljMLScaled, ljMLScaled, wl.first, wl.second));
    // find minimum and its index
    auto funcMinimum = std::min_element(funcValue.begin(), funcValue.end());
    MTindex minIndex = std::distance(funcValue.begin(), funcMinimum);
    
    if (minIndex < nRoots) {
        return roots[minIndex] / fUnitConv;
    } else {
        return lj_ML;
    }
}

MTfloat MutoMAP::surrogateFunctionValue(MTindex Mj, MTfloat ljn, MTfloat ljML, MTfloat lj, std::vector<MTfloat> w, std::vector<MTfloat> lm) {
    MTfloat value = 0.0;
    if (lj < 0) {
        return __DBL_MAX__;
    }
    // ML surrogate function
    value += 0 < lj ? static_cast<MTfloat>(Mj) * (std::log(lj) + ljML / lj)  : 0.0;
    // pair-wise lp-norm regularization surrogate function
    for (int i=0; i<w.size(); ++i) {
        value += 0.5 * fBeta * fP * w[i] * 
                 std::pow(std::pow(ljn-lm[i], 2) + fDelta*fDelta, (fP-2.0)/2.0) *
                 std::pow(2.0*lj - ljn - lm[i], 2);
    }
    
    // add object prior surroage function by choice
    if (fAddObjectPrior) {
        value += 0.5 * fBeta * fP * lj * lj *
                 std::pow(std::pow(ljn, 2) + fDelta*fDelta, (fP-2.0)/2.0);
    }
    return value;    
}

Matrix<MTindex, 3, 1> MutoMAP::index2xyz(MTindex i){
    Matrix<MTindex, 3, 1> xyz;
    // column major of Eigen default matrix/tensor
    xyz(0) = i % (fGrid.nx * fGrid.ny) % fGrid.nx;
    xyz(1) = i % (fGrid.nx * fGrid.ny) / fGrid.nx;
    xyz(2) = i / (fGrid.nx * fGrid.ny);
    return xyz;
}

std::vector<MTfloat> MutoMAP::getSixNeighbours(MTindex j, const Image& img) {
    std::vector<MTfloat> l;
    auto xyz = index2xyz(j);
    MTindex xp = xyz(0) < fGrid.nx-1 ? xyz(0) + 1 : fGrid.nx-1; // mirror boundary extension
    MTindex yp = xyz(1) < fGrid.ny-1 ? xyz(1) + 1 : fGrid.ny-1;
    MTindex zp = xyz(2) < fGrid.nz-1 ? xyz(2) + 1 : fGrid.nz-1;
    MTindex xm = 0 < xyz(0) ? xyz(0) - 1 : 0;
    MTindex ym = 0 < xyz(1) ? xyz(1) - 1 : 0;
    MTindex zm = 0 < xyz(2) ? xyz(2) - 1 : 0;
    l.push_back(img(xp, xyz(1), xyz(2)));
    l.push_back(img(xm, xyz(1), xyz(2)));
    l.push_back(img(xyz(0), yp, xyz(2)));
    l.push_back(img(xyz(0), ym, xyz(2)));
    l.push_back(img(xyz(0), xyz(1), zp));
    l.push_back(img(xyz(0), xyz(1), zm));
    for (auto &v : l) { v *= fUnitConv; } // scale values of lambdas
    return l;
}
std::vector<MTfloat> MutoMAP::getAllNeighbours(MTindex j, const Image& img) {
    std::vector<MTfloat> l;
    auto xyz = index2xyz(j);
    MTindex xp = xyz(0) < fGrid.nx-1 ? xyz(0) + 1 : fGrid.nx-1; // mirror boundary extension
    MTindex yp = xyz(1) < fGrid.ny-1 ? xyz(1) + 1 : fGrid.ny-1;
    MTindex zp = xyz(2) < fGrid.nz-1 ? xyz(2) + 1 : fGrid.nz-1;
    MTindex xm = 0 < xyz(0) ? xyz(0) - 1 : 0;
    MTindex ym = 0 < xyz(1) ? xyz(1) - 1 : 0;
    MTindex zm = 0 < xyz(2) ? xyz(2) - 1 : 0;
    // closest 
    l.push_back(img(xp, xyz(1), xyz(2)));
    l.push_back(img(xm, xyz(1), xyz(2)));
    l.push_back(img(xyz(0), yp, xyz(2)));
    l.push_back(img(xyz(0), ym, xyz(2)));
    l.push_back(img(xyz(0), xyz(1), zp));
    l.push_back(img(xyz(0), xyz(1), zm));
    // on the edge
    l.push_back(img(xp, yp, xyz(2)));
    l.push_back(img(xm, yp, xyz(2)));
    l.push_back(img(xp, ym, xyz(2)));
    l.push_back(img(xm, ym, xyz(2)));
    l.push_back(img(xyz(0), yp, zp));
    l.push_back(img(xyz(0), ym, zp));
    l.push_back(img(xyz(0), yp, zm));
    l.push_back(img(xyz(0), ym, zm));
    l.push_back(img(xp, xyz(1), zp));
    l.push_back(img(xm, xyz(1), zp));
    l.push_back(img(xp, xyz(1), zm));
    l.push_back(img(xm, xyz(1), zm));
    // at the corner
    l.push_back(img(xp, yp, zp));
    l.push_back(img(xm, yp, zp));
    l.push_back(img(xp, ym, zp));
    l.push_back(img(xm, ym, zp));
    l.push_back(img(xp, yp, zm));
    l.push_back(img(xm, yp, zm));
    l.push_back(img(xp, ym, zm));
    l.push_back(img(xm, ym, zm));
    for (auto &v : l) { v *= fUnitConv; } // scale values of lambdas
    return l;
}

// copy from old implementation of MLSDEM
MTindex MutoMAP::cubicEqnSolver(MTfloat* roots, MTfloat p3, MTfloat p2, MTfloat p1, MTfloat p0, MTfloat tol) {
    // get parameters based on the sign of p3
    MTfloat a, b, c; // for equation: x^3 + ax^2 + bx^ + c = 0
    if (std::abs(p3) < tol) {
        std::cout << "Warning: possible divide by zero in Cubic equation solver. p3 = " << p3 << std::endl;
    }
    a = p2 / p3;
    b = p1 / p3;
    c = p0 / p3;

    // initial value
	MTfloat r = 1 + std::max(std::max(std::abs(a), std::abs(b)), std::abs(c));
	MTfloat xinfl = -a / 3;									
	MTfloat gxinfl = xinfl * xinfl * xinfl + a * xinfl * xinfl + b * xinfl + c;
	MTfloat x0 = r;

	x0 = 0 < gxinfl ? -r : r;

	// solve for real roots by iteration
	MTfloat xk = x0;
	MTfloat rel = tol * 2;
	while (rel > tol) {
		MTfloat xk2 = xk * xk;
		MTfloat gxk = xk2 * xk + a * xk2 + b * xk + c;			// evaluation of func.
		MTfloat gpxk = 3 * xk2 + 2 * a * xk + b;			    // 1st order derivation
        MTfloat xknew = xk - gxk / gpxk;						// iteration result
		rel = std::abs((xknew - xk) / xk);
		xk = xknew;
	}	
	MTfloat x1 = xk;

	// get the other 2 possible real roots
	MTfloat D = x1 + a;
	MTfloat E = D * x1 + b;
	MTfloat signb = 1;
	signb = 0 <= D ? 1 : -1;

	MTfloat diak = D * D - 4 * E;

	// if they are complex numbers
	if (diak < 0) {
		roots[0] = x1;
		return 1;
	}
	
	MTfloat q = -0.5 * (D + signb * std::sqrt(diak));
	MTfloat x2 = q;
	MTfloat x3 = x2;
	x3 = tol < std::abs(q) ? E / q : x2;

	roots[0] = x1;
	roots[1] = x2;
	roots[2] = x3;
	return 3;
}