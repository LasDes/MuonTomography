#include "MutoSLR.h"
#include "MutoSiddon.h"

MutoSLR::MutoSLR() {
    // default voxel world
    fGrid.nx = fGrid.ny = fGrid.nz = 100;
    fGrid.x_min = fGrid.y_min = fGrid.z_min = 0.0;
    fGrid.dx = fGrid.dy = fGrid.dz = 10.0;
    // default settings for add weight to scattering point
    fLengthWeighted = false;
    fHasEnergy = false;
    fLogPrint = false;
    // reset counts
    fNTotal = 0;
    // scattering angle unit, default to rad, i.e. x1000.0 for mrad
    fAngMag = 1.0;
    // eps for floating points comparison
    fEPS = 1.0e-8;
}

MutoSLR::MutoSLR(json config) : MutoSLR() {
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
    fLengthWeighted = config.value("length_weighted", false);
    fHasEnergy = config.value("has_energy", false);
    fLogPrint = config.value("log_results", false);
    fNTotal = config.value("total_num_rays", 0);
    fAngMag = config.value("angle_unit_mag", 1.0);
    fEPS = config.value("eps", 1.0e-8);
}

MutoSLR::~MutoSLR() {}

Image MutoSLR::reconstruct(const MutoMuonData& rays) {
    Image img = initailizeImage(); // contain scat. angles
    Image imgWeight = initailizeImage(); // contain weight of voxels
    MutoSiddon siddon(fGrid);

    fNTotal = 0 < fNTotal ? fNTotal : rays.size();

    // setting minimum and maximum values for height of voxel region
    MTfloat zVoxelMin = fGrid.z_min;
    MTfloat zVoxelMax = fGrid.z_min + fGrid.dz * static_cast<MTfloat>(fGrid.nz);

    // loop each muon in the data
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
        
        for (auto v : voxels) {
            MTfloat iValue = 0.0;
            if (fLengthWeighted) {
                imgWeight(v.x, v.y, v.z) += v.length / fGrid.dz;
                iValue = (axz*axz + ayz*ayz) * 0.5 * (v.length/lTotal)*(v.length/lTotal); 
            } else {
                imgWeight(v.x, v.y, v.z) += 1.0;
                iValue = (axz*axz + ayz*ayz) * 0.5 / static_cast<MTfloat>(voxels.size() * voxels.size());
            }
            if (fHasEnergy) {
                iValue *= rays[i].ein * rays[i].ein;
            }
            img(v.x, v.y, v.z) += iValue;
        } // end of adding weight to passing voxels
    } // end of muon rays loop

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

Image MutoSLR::initailizeImage() {
    Image init(fGrid.nx, fGrid.ny, fGrid.nz);
    init.setZero(); // set all elements to zeros
    return init;
}

void MutoSLR::logResults() {
    std::cout << "==================================================" << std::endl;
    std::cout << "|| SLR reconstruction results                  ||" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "number of rays processed: " << fNTotal << std::endl;
}
