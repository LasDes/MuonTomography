#include "MutoPoCA.h"
#include "MutoSiddon.h"

MutoPoCA::MutoPoCA() {
    // default voxel world
    fGrid.nx = fGrid.ny = fGrid.nz = 100;
    fGrid.x_min = fGrid.y_min = fGrid.z_min = 0.0;
    fGrid.dx = fGrid.dy = fGrid.dz = 10.0;
    // default settings for add weight to scattering point
    fLengthWeighted = true;
    fStraightLine = false;
    fHasEnergy = false;
    fLogPrint = false;
    // reset counts
    fNTotal = 0;
    fNNotParallel = 0;
    fNScatInGrid = 0;
    // scattering angle unit, default to mrad, i.e. x1000.0
    fAngMag = 1.0;
    // eps for floating points comparison
    fEPS = 1.0e-8;
}

MutoPoCA::MutoPoCA(json config) : MutoPoCA() {
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
    fLengthWeighted = config.value("length_weighted", true);
    fStraightLine = config.value("straight_line", false);
    fHasEnergy = config.value("has_energy", false);
    fLogPrint = config.value("log_results", false);
    fNTotal = config.value("total_num_rays", 0);
    fAngMag = config.value("angle_unit_mag", 1.0);
    fEPS = config.value("eps", 1.0e-8);
}

MutoPoCA::~MutoPoCA() {}

Image MutoPoCA::reconstruct(const MutoMuonData& rays) {
    Image img = initailizeImage(); // contain scat. angles
    Image imgWeight = initailizeImage(); // contain weight of voxels
    MutoSiddon siddon(fGrid);

    fNTotal = 0 < fNTotal ? fNTotal : rays.size();
    fNNotParallel = 0; // clear counters
    fNScatInGrid = 0;

    // setting minimum and maximum values for z (scattering outside of detection area)
    MTfloat zDetMin = rays[0].pout(2);
    MTfloat zDetMax = rays[0].pin(2);

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
        // add scattering angle to voxel that contains pScat
        if (fHasEnergy) { // not implemented
            img(xscat, yscat, zscat) += aScatSquare * rays[i].ein * rays[i].ein;
        } else {
            img(xscat, yscat, zscat) += aScatSquare;
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
        logResults();
    }

    // processs image with weight matrix and return
    for (MTindex i=0; i<img.size(); ++i) {
        img(i) = imgWeight(i) < fEPS ? 0 : img(i) / imgWeight(i);
    }
    return img;
}

Image MutoPoCA::initailizeImage() {
    Image init(fGrid.nx, fGrid.ny, fGrid.nz);
    init.setZero(); // set all elements to zeros
    return init;
}

Vector3 MutoPoCA::getPointOfClosestApproach(const RayData& ray) {
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
    fNNotParallel ++;
    // calculate coefficients in both lines
    Cin = (b*e - c*d) / denominant;
    Cout = (a*e - b*d) / denominant;

    return 0.5 * (ray.pin + Cin*ray.din + ray.pout + Cout*ray.dout);    
}

MTfloat MutoPoCA::getScatAngleSquare(const Vector3& d1, const Vector3& d2) {
    // NOTE: an alternate way to calcuate scattering angle is to calculate 
    //       the dot product of two directions (cos angle)
    // MTfloat cosA = d1.dot(d2) / d1.norm() / d2.norm();
    // return std::pow(std::acos(cosA) * fAngMag, 2);

    MTfloat axz = (std::atan2(d1(0), d1(2)) - std::atan2(d2(0), d2(2))) * fAngMag;
    MTfloat ayz = (std::atan2(d1(1), d1(2)) - std::atan2(d2(1), d2(2))) * fAngMag;
    return (axz * axz + ayz * ayz) / 2.0; // return average value of ax^2 and ay^2
}

void MutoPoCA::logResults() {
    std::cout << "==================================================" << std::endl;
    std::cout << "|| PoCA reconstruction results                  ||" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "number of rays processed: " << fNTotal << std::endl;
    std::cout << "number of unparalleled rays: " << fNNotParallel << std::endl;
    std::cout << "number of interactions in voxels: " << fNScatInGrid << std::endl;
}
