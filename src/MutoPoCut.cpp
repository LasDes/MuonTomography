#include "MutoPoCut.h"
#include "MutoSiddon.h"

MutoPoCut::MutoPoCut() : MutoPoCA() {
    // fix some parameters
    fStraightLine = false;
    
    // new parameter
    fPercentile = 1.0;
}

MutoPoCut::MutoPoCut(json config) : MutoPoCA(config) {
    fStraightLine = false;
    fPercentile = config.value("cut_percentile", 1.0);
}

MutoPoCut::~MutoPoCut() {}

Image MutoPoCut::reconstruct(const MutoMuonData& rays) {
    std::cout << "==================================================" << std::endl;
    std::cout << "|| PoCA with threshold for scattering angles    ||" << std::endl;
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
        MTfloat aScatSquare = getScatAngleSquare(rays[i].din, rays[i].dout);
        
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
        // add scattering angle to voxel's scattering store
        if (fHasEnergy) { // not implemented
            imgScat[xyz2index(xscat, yscat, zscat)].push_back(aScatSquare * rays[i].ein * rays[i].ein);
        } else {
            imgScat[xyz2index(xscat, yscat, zscat)].push_back(aScatSquare);
        }  

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
        std::cout << "Start assigning scattering density to each voxel with threshold. " << std::endl;
    }

    // processs image with weight matrix and return
    for (MTindex i=0; i<img.size(); ++i) {
        std::sort(imgScat[i].begin(), imgScat[i].end());
        MTindex threshold = static_cast<MTindex>(imgScat[i].size() * fPercentile);
        MTfloat imgValue = std::accumulate(imgScat[i].begin(), imgScat[i].begin()+threshold, 0.0);
        img(i) = imgWeight(i) < fEPS ? 0 : imgValue / imgWeight(i);
    }

    // log if desired
    if (fLogPrint) {
        logResults();
    }

    return img;
}
