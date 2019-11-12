#include "MutoSLRIter.h"
#include "MutoSiddon.h"

MutoSLRIter::MutoSLRIter() : MutoSLR() {
    // default number of iteration
    fNIteration = 1;
}

MutoSLRIter::MutoSLRIter(json config) : MutoSLR(config) {
    fNIteration = config.value("num_iteration", 1);
}

MutoSLRIter::~MutoSLRIter() {}

Image MutoSLRIter::reconstruct(const MutoMuonData& rays) {
    Image img = initailizeImage(); // contain scat. angles
    Image imgWeight = initailizeImage(); // contain weight of voxels
    MutoSiddon siddon(fGrid);

    fNTotal = 0 < fNTotal ? fNTotal : rays.size();

    // setting minimum and maximum values for height of voxel region
    MTfloat zVoxelMin = fGrid.z_min;
    MTfloat zVoxelMax = fGrid.z_min + fGrid.dz * static_cast<MTfloat>(fGrid.nz);

    // loop each muon in the data
    for (int iter=0; iter<fNIteration; ++iter) {
        for (int i=0; i<fNTotal; ++i) {
            // calculate scattering information from current muon
            MTfloat axz = (std::atan2(rays[i].din(0), rays[i].din(2)) - std::atan2(rays[i].dout(0), rays[i].dout(2))) * fAngMag;
            MTfloat ayz = (std::atan2(rays[i].din(1), rays[i].din(2)) - std::atan2(rays[i].dout(1), rays[i].dout(2))) * fAngMag;

            // use a straight line to approximate muon path
            MTfloat tIn = (zVoxelMax - rays[i].pin(2)) / rays[i].din(2);
            MTfloat tOut = (zVoxelMin - rays[i].pout(2)) / rays[i].dout(2);
            Vector3 pIncoming = rays[i].pin + tIn * rays[i].din;
            Vector3 pOutgoing = rays[i].pout + tOut * rays[i].dout;

            // std::cout << "Inray: " << rays[i].pin << " | " << rays[i].din << std::endl;
            // std::cout << "Outray: " << rays[i].pout << " | " << rays[i].dout << std::endl;
            // std::cout << pIncoming << " | " << pOutgoing << std::endl;

            auto voxels = siddon.getVoxelPath(pIncoming, pOutgoing);
            MTfloat lTotal = (pIncoming - pOutgoing).norm();
            MTfloat valTotal = 0.0;
            for (auto v : voxels) {
                valTotal += img(v.x, v.y, v.z);
            }
            
            for (auto v : voxels) {
                MTfloat iValue = 0.0;
                if (fLengthWeighted) {
                    imgWeight(v.x, v.y, v.z) += v.length / fGrid.dz;
                    iValue = (axz*axz + ayz*ayz) * 0.5 * (v.length/lTotal)*(v.length/lTotal); 
                } else {
                    imgWeight(v.x, v.y, v.z) += 1.0;
                    iValue = (axz*axz + ayz*ayz) * 0.5 / static_cast<MTfloat>(voxels.size() * voxels.size());
                }
                // weight scattering angles by current density values
                if (iter != 0) {
                    iValue *= img(v.x, v.y, v.z) / valTotal;
                }

                if (fHasEnergy) {
                    iValue *= rays[i].ein * rays[i].ein;
                }
                img(v.x, v.y, v.z) += iValue;
            } // end of adding weight to passing voxels
        } // end of muon rays loop
        std::cout << "Iteration #" + std::to_string(iter) + " finished. " << std::endl;
    } // end of iteration loop

    // log if desired
    if (fLogPrint) {
        logResults();
    }

    // processs image with weight matrix and return
    for (MTindex i=0; i<img.size(); ++i) {
        img(i) = imgWeight(i) < fEPS ? 0 : img(i) / imgWeight(i);
    }
    return img;
}
