#include "MutoMLSD.h"
#include "MutoFile.h"
#include "MutoSiddon.h"
#include <random>
// #include <cassert>

// --------------------------------------------------------------------
// Constructors and destructor
// --------------------------------------------------------------------
MutoMLSD::MutoMLSD(){
    // default voxel world
    fGrid.nx = fGrid.ny = fGrid.nz = 100;
    fGrid.x_min = fGrid.y_min = fGrid.z_min = 0.0;
    fGrid.dx = fGrid.dy = fGrid.dz = 10.0;

    // default settings for add weight to scattering point
    fUpdateWithMean = false;
    fStraightLine = false;
    fHasEnergy = false;
    fEarlyStop = false;
    fLogPrint = false;

    // reset counts
    fNTotal = 0;
    fNScatInGrid = 0;

    // mlsd iteration parameters
    fPercentileMean = 100;
    fNIteration = 10;

    // initialize methods
    fInitMethod = 0; // 0: image, 1: random, 2: fixed value
    fInitValue = 1.0e-7;
    fInitImageFile = ""; // default empty path for initial image

    // error matrix estimation for real measurement system
    fErrorXY = 1.0;
    fErrorZo = 1000.0;
    fErrorZi = 1000.0;

    // scattering angle unit, default to mrad, i.e. x1000.0
    fAngMag = 1.0;
    // eps for floating points comparison
    fEPS = 1.0e-8;
}

MutoMLSD::MutoMLSD(json config) : MutoMLSD() {
    // use json `value(key, default)` to set internal parameters
    try {
        json grid = config.at("grid");
        fGrid.nx = grid.value("nx", 100);
        fGrid.ny = grid.value("ny", 100);
        fGrid.nz = grid.value("nz", 100);
        fGrid.x_min = grid.value("x_min", 0.0);
        fGrid.y_min = grid.value("y_min", 0.0);
        fGrid.z_min = grid.value("z_min", 0.0);
        fGrid.dx = grid.value("dx", 10.0);
        fGrid.dy = grid.value("dy", 10.0);
        fGrid.dz = grid.value("dz", 10.0);
    } catch(json::out_of_range) {}

    // other setting are in the some level as grid
    fUpdateWithMean = config.value("update_with_mean", false);
    fStraightLine = config.value("straight_line", false);
    fHasEnergy = config.value("has_energy", false);
    fEarlyStop = config.value("early_stop", false);
    fLogPrint = config.value("log_results", false);
    fNTotal = config.value("total_num_rays", 0);
    fPercentileMean = config.value("mean_percentile", 100);
    fNIteration = config.value("num_iteration", 10);

    // initializations
    fInitMethod = config.value("initialize_method", 0);
    fInitValue = config.value("initialize_value", 1.0e-7);
    fInitImageFile = config.value("initial_image_file", "");

    // error matrix 
    fErrorXY = config.value("error_xy", 1.0);
    fErrorZo = config.value("error_dz_outer", 1000.0);
    fErrorZi = config.value("error_dz_inner", 1000.0);

    fAngMag = config.value("angle_unit_mag", 1.0);
    fEPS = config.value("eps", 1.0e-8);
}

MutoMLSD::~MutoMLSD() {}

// --------------------------------------------------------------------
// Core interface of MLSD method
// --------------------------------------------------------------------
Image MutoMLSD::reconstruct(const MutoMuonData& rays) {
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

// --------------------------------------------------------------------
// Private functions for initialization
// --------------------------------------------------------------------

Image MutoMLSD::initailizeImage() {
    // generate random noise as initial image
    if (fInitMethod == 0) {
        auto imgAndHeader = MutoFile::loadImage(fInitImageFile);
        auto inGrid = imgAndHeader.second.grid;
        if (imgAndHeader.first.size() == 0) {
            std::cout << "Warning: Initial image file <" << fInitImageFile << "> not exist. " << std::endl;
            std::cout << "\tSetting initialization method to random numbers ..." << std::endl;
            fInitMethod = 1;
        } else if (inGrid.nx != fGrid.nx ||  inGrid.ny != fGrid.ny || inGrid.nz != fGrid.nz ) { // should implement equality for VoxelGrid, simplified here.
            // future feature: scale the image to desired dimensions
            // use input image's grid 
            std::cout << "Warning: Initial image's dimensions do not match the configuration. " << std::endl;
            std::cout << "\tResetting grid configuration as defined in input image  ..." << std::endl;
            fGrid = inGrid;
            return imgAndHeader.first;
        } else {
            return imgAndHeader.first;
        }
    }

    // other initlization methods
    Image img (fGrid.nx, fGrid.ny, fGrid.nz);

    // random number generators
    std::random_device rdev;
    std::default_random_engine generator{rdev()};
    std::normal_distribution<double> distribution(2.0, 0.5); // default to value close to water

    if (fInitMethod == 1) {
        // TODO: change to TensorBase -> setRandom to improve efficiency
        MTfloat lTemp;
        for (MTindex i=0; i<img.size(); ++i) {
            lTemp = distribution(generator) * 1.0e-7;
            img(i) = 1.0e-10 <= lTemp ? lTemp : 1.0e-10;
        }
    } else { // set image to a fixed value
        img.setConstant(fInitValue);
    }
    return img;
}

void MutoMLSD::gatherInformation(const MutoMuonData & data, std::vector<MTindex>& index, 
                                 std::vector<MTfloat>& Sx, std::vector<MTfloat>& Dx, 
                                 std::vector<MTfloat>& Sy, std::vector<MTfloat>& Dy, 
                                 std::vector<MTfloat>& p, std::vector<Vector3>& pScatVec){
    MTindex iMuon = 0;
    MTfloat zDetMin = data[0].pout(2);
    MTfloat zDetMax = data[0].pin(2);
    MTfloat pNominal = 3.0; // in the unit of GeV

    // reserve memory assuming that half of the muons scattered inside object
    index.reserve(data.size() / 2);
    Sx.reserve(data.size() / 2);
    Dx.reserve(data.size() / 2);
    Sy.reserve(data.size() / 2);
    Dy.reserve(data.size() / 2);
    pScatVec.reserve(data.size() / 2);

    if (fHasEnergy) {
        p.reserve(data.size() / 2);
    } else { // assign a constant value to p
        p.push_back(pNominal / pNominal); // in old MLSDEM implementation, this value is 3/2 for 2GeV muon
    }

    for (auto ray : data) {
        // confirm index of current muon with preset fNTotal
        if (iMuon == fNTotal) { break; }

        Vector3 pScat = getPointOfClosestApproach(ray);

        // check scattering in the grid system and detection area
        if (pScat(0) < fGrid.x_min || fGrid.x_min + fGrid.dx * fGrid.nx <= pScat(0) || 
            pScat(1) < fGrid.y_min || fGrid.y_min + fGrid.dy * fGrid.ny <= pScat(1) || 
            pScat(2) < fGrid.z_min || fGrid.z_min + fGrid.dz * fGrid.nz <= pScat(2) ||
            pScat(2) < zDetMin || zDetMax < pScat(2)) { // out of detection zone
            iMuon ++;
            continue;
        }

        // push scattering points
        index.push_back(iMuon);
        pScatVec.push_back(pScat);
        fNScatInGrid ++;

        // get scattering angles and displacements in XZ and YZ plane resp.
        auto vscat = getScatteringAngle(ray.din, ray.dout);
        auto vdisp = getDisplacement(ray);
        Sx.push_back(vscat.first);
        Sy.push_back(vscat.second);
        Dx.push_back(vdisp.first);
        Dy.push_back(vdisp.second);

        // calculate momentum if energy is available
        if (fHasEnergy) {
            p.push_back(pNominal / ray.ein); // definition: pr = p0 / p
        }

        iMuon ++;
    }
}

Vector3 MutoMLSD::getPointOfClosestApproach(const RayData& ray) {
    // based on Online resource about distance btw lines: http://geomalgorithms.com/a07-_distance.html
    MTfloat a, b, c, d, e; // intermediate values
    MTfloat Cin, Cout; // coefficient of closest point in both lines
    Vector3 w0 = ray.pin - ray.pout; // vector btw two reference points
    a = ray.din.dot(ray.din);
    b = ray.din.dot(ray.dout);
    c = ray.dout.dot(ray.dout);
    d = ray.din.dot(w0);
    e = ray.dout.dot(w0);
    // check for paralleled lines
    MTfloat denominant = a*c - b*b;
    if (std::abs(denominant) < fEPS) {
        return 0.5 * (ray.pin + ray.pout);
    }
    // calculate coefficients in both lines
    Cin = (b*e - c*d) / denominant;
    Cout = (a*e - b*d) / denominant;

    return 0.5 * (ray.pin + Cin*ray.din + ray.pout + Cout*ray.dout);    
}

std::pair<MTfloat, MTfloat> MutoMLSD::getScatteringAngle(const Vector3& d1, const Vector3& d2) {
    MTfloat axz = (std::atan2(d1(0), d1(2)) - std::atan2(d2(0), d2(2))) ; //* fAngMag;
    MTfloat ayz = (std::atan2(d1(1), d1(2)) - std::atan2(d2(1), d2(2))) ; //* fAngMag;
    return std::pair<MTfloat, MTfloat> (axz, ayz);
}

std::pair<MTfloat, MTfloat> MutoMLSD::getDisplacement(const RayData& ray) {
    // in and out angle in XZ and YZ plane (0:x, 1:y, 2:z)
    MTfloat axin = std::atan2(ray.din(0), ray.din(2));
    MTfloat axout = std::atan2(ray.dout(0), ray.dout(2));
    MTfloat ayin = std::atan2(ray.din(1), ray.din(2));
    MTfloat ayout = std::atan2(ray.dout(1), ray.dout(2));
    // project incoming / outgoing ray to z_min of the grid system. (p: projected, t: true)
    Vector3 pProj = ray.pin + (fGrid.z_min - ray.pin(2)) / ray.din(2) * ray.din;
    Vector3 pTrue = ray.pout + (fGrid.z_min - ray.pout(2)) / ray.dout(2) * ray.dout;
    MTfloat xp = pProj(0);
    MTfloat xt = pTrue(0);
    MTfloat yp = pProj(1);
    MTfloat yt = pTrue(1);

    // calculate displacement 
    MTfloat Lxy = std::sqrt(1 + std::pow(std::tan(axin), 2) + std::pow(std::tan(ayin), 2));
    MTfloat dxz = (xt - xp) * std::cos(axin) * Lxy * std::cos(axout) / std::cos(axout - axin);
    MTfloat dyz = (yt - yp) * std::cos(ayin) * Lxy * std::cos(ayout) / std::cos(ayout - ayin);

    return std::pair<MTfloat, MTfloat> (dxz, dyz);
}

// --------------------------------------------------------------------
// Private functions for MLSD iteration
// --------------------------------------------------------------------

MTindex MutoMLSD::getMuonPathAndTail(const RayData& ray, const Vector3& pScat, std::vector<VoxelData>& path, 
                                     std::vector<MTfloat>& tail, MutoSiddon& siddon) {
    // clear old data
    if (path.size() != 0) { path.clear(); } // check size before clear 
    if (tail.size() != 0) { tail.clear(); } // this line is not necessary, tail vector will be resized and loop over

    // use siddon algorithm to get all the passed voxels
    if (fStraightLine) {
        auto voxels = siddon.getVoxelPath(ray.pin, ray.pout);
        path.insert(path.begin(), voxels.begin(), voxels.end());
    } else {
        auto voxelsIn = siddon.getVoxelPath(ray.pin, pScat);
        auto voxelsOut = siddon.getVoxelPath(pScat, ray.pout);
        // take care of the scattering voxel
        voxelsIn.back().length += voxelsOut[0].length;
        path.insert(path.begin(), voxelsIn.begin(), voxelsIn.end());
        path.insert(path.end(), voxelsOut.begin()+1, voxelsOut.end());
        // assert(path.size() == voxelsIn.size() + voxelsOut.size() -1);
    }

    // get tailing length of each voxel
    tail.resize(path.size());
    tail.back() = 0.0;
    if (1 < tail.size()) { // loop backwards to get tailing lengths
        for (int ri=tail.size()-2; 0<=ri; --ri) {
            tail[ri] = tail[ri+1] + path[ri+1].length;
        }
    }
    return path.size();
}

MTindex MutoMLSD::generateWeightMatrix(std::vector<std::pair<MTindex, Mat2x2>>& W, 
                                       const std::vector<VoxelData>& v, const std::vector<MTfloat>& t){
    // unit weight matrix
    Mat2x2 Wij;
    MTfloat LT;
    for (MTindex j=0; j<t.size(); ++j) {
        LT = v[j].length * t[j];
        Wij(0,0) = v[j].length;
        Wij(0,1) = Wij(1,0) = LT + v[j].length * v[j].length * 0.5;
        Wij(1,1) = v[j].length * v[j].length * v[j].length / 3.0 + (v[j].length + t[j]) * LT;
        W.push_back(std::pair<MTindex, Mat2x2>(xyz2index(v[j].x, v[j].y, v[j].z), Wij));
    }
    return W.size();
}

void MutoMLSD::allocateSValueMatrix(std::vector<std::vector<MTfloat>>& S, const std::vector<VoxelData>& path){
    for (auto v : path) {
        S[xyz2index(v.x, v.y, v.z)].push_back(0.0);
    }
}

// --------------------------------------------------------------------
// Private functions for MLSD iteration
// --------------------------------------------------------------------

Mat2x2 MutoMLSD::getCovarianceMatrix(const Image& img, const std::vector<std::pair<MTindex, Mat2x2>>& W, MTfloat p){
    
    Mat2x2 cov;
    cov.setZero();
    // add values to cov
    for (auto Wj : W) { cov += img(Wj.first) * Wj.second; }
    return p * p * cov;
}

MTfloat MutoMLSD::updateSValueMatrix(std::vector<std::vector<MTfloat>>& S, std::vector<MTindex>& Sindex, const Image&img, 
                                     const Mat2x2& cov, const std::vector<std::pair<MTindex, Mat2x2>>& W, 
                                     MTfloat p, MTfloat Sx, MTfloat Sy, MTfloat Dx, MTfloat Dy) {
    // calculate inverse of covariance matrix
    MTfloat cost = 0.0;
    MTfloat SigmaDet = 1.0;
    Mat2x2 SigmaInv;
    SigmaInv.setZero();
    bool covInversible;
    cov.computeInverseAndDetWithCheck(SigmaInv, SigmaDet, covInversible);
    if (!covInversible || SigmaDet <= 0) {
        std::cout << "Warning: one Matrix is not inversible. " << std::endl;
        cost = 0.0; // prevent divide by zero error
    } else { // calculate cost 
        cost = std::log(SigmaDet) + 
               ((Sx*Sx + Sy*Sy) * SigmaInv(0,0) + 
                (Dx*Dx + Dy*Dy) * SigmaInv(1,1) + 
                (Sx*Dx + Sy*Dy) * 2.0 * SigmaInv(0,1)) * 0.5;
    }
    // update S matrix for all voxel j that passed by Muon Ray i
    for (auto Wj : W) {
        MTindex j = Wj.first;
        Mat2x2 wij = Wj.second;
        MTfloat lj = img(j);
        MTfloat pr2lj2 = p * p * lj * lj;
        MTfloat xySij;
        MTfloat trCovW = (SigmaInv * wij).trace();
        Mat2x2 CovWCov = SigmaInv * wij * SigmaInv;

        xySij = (Sx*Sx + Sy*Sy) * CovWCov(0,0) + 
                (Dx*Dx + Dy*Dy) * CovWCov(1,1) + 
                (Sx*Dx + Sy*Dy) * (CovWCov(0,1) + CovWCov(1,0));
        MTfloat sij = 2.0 * lj + (0.5 * xySij - trCovW) * pr2lj2;

        // update S_values (zero values included based on the original paper)
        S[j][Sindex[j]] = sij;
        Sindex[j] ++;  
    }
    // return cost for this muon
    return cost;
}

MTfloat MutoMLSD::updateImageVoxel(std::vector<MTfloat>& Sj, MTindex Mj) {
    if (Mj == 0) {
        return 0.0;
    }
    if (fUpdateWithMean) {
        if (fPercentileMean == 100) {
            return std::accumulate(Sj.begin(), Sj.begin()+Mj, 0.0) / static_cast<MTfloat>(Mj) * 0.5;
        } else {
            MTindex jMean = static_cast<MTindex>(fPercentileMean/100.0*Mj);
            if (jMean < 1) jMean = 1; // prevent zero values in division
            std::sort(Sj.begin(), Sj.end());
            return std::accumulate(Sj.begin(), Sj.begin()+jMean, 0.0) / static_cast<MTfloat>(jMean) * 0.5;
        }
    } else { // update using median value of Sj
        // two cases: Mj is even or odd number
        std::nth_element(Sj.begin(), Sj.begin() + Mj/2, Sj.begin() + Mj);
        MTfloat median = Sj[Mj/2];
        if (Mj % 2 == 0) { // find [Mj/2] and [Mj/2+1]
            std::nth_element(Sj.begin(), Sj.begin() + Mj/2-1, Sj.begin() + Mj);
            return (median + Sj[Mj/2-1]) * 0.25;
        } else {
            return median * 0.5;
        }        
    }
}


// --------------------------------------------------------------------
// Private functions for print information
// --------------------------------------------------------------------

void MutoMLSD::logMessage(const std::string& msg, bool printTime) {
    if (fLogPrint) {
        std::cout << msg;
        if (printTime) {
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - fTimeStart);
            std::cout <<  " Time: " << time_span.count() << " seconds. ";
        }
        std::cout << std::endl;
    }
}

void MutoMLSD::logResults() {
    std::cout << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "|| MLSD reconstruction results                  ||" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "number of rays processed: " << fNTotal << std::endl;
    std::cout << "number of interactions in voxels: " << fNScatInGrid << std::endl;
    std::cout << "number of iterations: " << fNIteration << std::endl;
}
