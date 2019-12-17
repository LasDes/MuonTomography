#include "MutoPoCAGMM.h"
#include "MutoSiddon.h"

MutoPoCAGMM::MutoPoCAGMM() : MutoPoCA() {
    // initial values for some parameters
    const double k[] = {45.07194375542501,
                        13.343386349742266,
                        6.336775578006414,
                        3.480783662355906,
                        2.0086482672843475,
                        1.1684165786047735,
                        0.6633338267301239,
                        0.34558879777808715,
                        0.15154595131997278,
                        0.035501855168395866};
    fKList.insert(fKList.end(), std::begin(k), std::end(k));
    fTolerance = 1.0e-4; // relative tolerance for update 
    fMaxIteration = 100; // max iteration steps
}

MutoPoCAGMM::MutoPoCAGMM(json config) : MutoPoCA(config) {
    fTolerance = config.value("ML_tolerance", 1.0e-4);
    fMaxIteration = config.value("ML_maxiter", 100);
    try {
        auto k = config.at("ML_multiplier");
        fKList.clear();
        fKList.insert(fKList.end(), std::begin(k), std::end(k));
    } catch(json::out_of_range) {}
}

MutoPoCAGMM::~MutoPoCAGMM() {}

Image MutoPoCAGMM::reconstruct(const MutoMuonData& rays) {
    std::cout << "==================================================" << std::endl;
    std::cout << "|| PoCA with Gaussian Mixture Model Fit         ||" << std::endl;
    std::cout << "==================================================" << std::endl;
    Image img = initailizeImage(); // contain scat. angles
    Image imgWeight = initailizeImage(); // contain weight of voxels
    MutoSiddon siddon(fGrid);

    // scattering angle store
    std::vector<std::vector<MTfloat>> imgScat;
    imgScat.resize(img.size());

    fNTotal = 0 < fNTotal ? fNTotal : rays.size();
    fNNotParallel = 0; // clear counters
    fNScatInGrid = 0;

    // setting minimum and maximum values for z (scattering outside of detection area)
    MTfloat zDetMin = rays[0].pout(2);
    MTfloat zDetMax = rays[0].pin(2);

    auto xyz2index = [&](MTindex x, MTindex y, MTindex z) { return x + y*fGrid.nx + z*fGrid.nx*fGrid.ny; };

    // loop each muon in the data
    for (int i=0; i<fNTotal; ++i) {
        // calculate scattering information from current muon
        Vector3 pScat = getPointOfClosestApproach(rays[i]);
        MTfloat axz = (std::atan2(rays[i].din(0), rays[i].din(2)) - std::atan2(rays[i].dout(0), rays[i].dout(2))) * fAngMag;
        MTfloat ayz = (std::atan2(rays[i].din(1), rays[i].din(2)) - std::atan2(rays[i].dout(1), rays[i].dout(2))) * fAngMag;
        
        // get index of pScat
        MTindex xscat = static_cast<MTindex>((pScat(0) - fGrid.x_min) / fGrid.dx);
        MTindex yscat = static_cast<MTindex>((pScat(1) - fGrid.y_min) / fGrid.dy);
        MTindex zscat = static_cast<MTindex>((pScat(2) - fGrid.z_min) / fGrid.dz);
        // check scattering in the voxel world
        if (xscat < 0 || fGrid.nx <= xscat || 
            yscat < 0 || fGrid.ny <= yscat || 
            zscat < 0 || fGrid.nz <= zscat ||
            pScat(2) < zDetMin || zDetMax < pScat(2)) { // out of detection zone
            continue;
        } else {
            fNScatInGrid ++;
        }
        // add scattering angle (in xz and yz planes) to voxel's scattering store
        imgScat[xyz2index(xscat, yscat, zscat)].push_back(axz);
        imgScat[xyz2index(xscat, yscat, zscat)].push_back(ayz);

        // use a straight line to approximate muon path
        if (fStraightLine) {
            auto voxels = siddon.getVoxelPath(rays[i].pin, rays[i].pout);
            for (auto v : voxels) {
                if (fLengthWeighted) {
                    imgWeight(v.x, v.y, v.z) += v.length;
                } else {
                    imgWeight(v.x, v.y, v.z) += fGrid.dz;
                }
            }
        } else { // use two parts of muon path centered at pScat
            // incoming and outgoing rays
            auto voxelsIn = siddon.getVoxelPath(rays[i].pin, pScat);
            auto voxelsOut = siddon.getVoxelPath(pScat, rays[i].pout);

            for (auto v : voxelsIn) {
                if (fLengthWeighted) {
                    imgWeight(v.x, v.y, v.z) += v.length;
                } else {
                    imgWeight(v.x, v.y, v.z) += fGrid.dz;
                }
            }
            for (auto v : voxelsOut) {
                if (fLengthWeighted) {
                    imgWeight(v.x, v.y, v.z) += v.length;
                } else {
                    if (v.x != xscat || v.y != yscat || v.z != zscat){
                        imgWeight(v.x, v.y, v.z) += fGrid.dz;
                    }  // prevent doubling the weight at the scattering voxel
                }
            }
        } // end of adding weight to passing voxels
    } // end of muon rays loop

    // log if desired
    if (fLogPrint) {
        std::cout << "Finish processing each muon. " << std::endl;
        std::cout << "Start assigning scattering density to each voxel by fitting to GMM. " << std::endl;
    }

    // processs image with weight matrix and return
    for (MTindex i=0; i<img.size(); ++i) {
        MTindex dataSize = imgScat[i].size();
        MTfloat varValue = std::accumulate(imgScat[i].begin(), imgScat[i].end(), 0.0, 
                                           [](MTfloat s, MTfloat v){return s + v*v;});

        MTfloat linit = imgWeight(i) < fEPS ? 0 : varValue / imgWeight(i) / (9/1.075*1.075);
        img(i) = fitGMM(imgScat[i], std::max(fEPS, std::min(linit, 1.0e-4)));
    }

    // log if desired
    if (fLogPrint) {
        logResults();
    }

    return img;
}

MTfloat MutoPoCAGMM::fitGMM(std::vector<MTfloat> data, MTfloat linit) {
    return 0;
}

/*
def gaussian_prob(x, theta):
    theta_k = (k * theta).reshape(1,-1)
    x = x.reshape(-1,1)
    w = 1.0 / np.sqrt(2*np.pi*theta_k) * np.exp(-x**2 / theta_k / 2.0)
    sw = w.sum(axis=1).reshape(-1,1)
    sw[sw < eps] = eps
    return w / sw
    
def update_theta(x, w):
    nk = w.sum(axis=0).reshape(1,-1)
    return (np.dot(x.reshape(1,-1)**2, w) / nk / k).mean()

def EMSolver(x, theta0, rlimit=1e-4, maxiter=100):
    l_last = theta0
    w = gaussian_prob(x, l_last)
    for i in range(maxiter):
        l = update_theta(x, w)
        if np.abs(l - l_last) / l_last < rlimit:
            break
        l_last = l
        w = gaussian_prob(x, l)
    return l, i
*/