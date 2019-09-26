#include "MutoSiddon.h"
#include <cmath>

MutoSiddon::MutoSiddon(const VoxelGrid & grid) {
    fGrid = grid;
}

std::vector<VoxelData> MutoSiddon::getVoxelPath(const Vector3 & p1, const Vector3 & p2) {
    // initialize return data
    std::vector<VoxelData> vpath;

    // get first and last intersection points for each parallel planes
    // i.e. amin and amax for x, y, z 
    MTfloat a0x = 0.0, a0y = 0.0, a0z = 0.0, anx = 1.0, any = 1.0, anz = 1.0;
    MTfloat eps = 1.0e-4;

    if (std::abs(p2(0) - p1(0)) > eps) {
        a0x = (fGrid.x_min - p1(0)) / (p2(0) - p1(0));
        anx = (fGrid.x_min + fGrid.dx * static_cast<MTfloat>(fGrid.nx) - p1(0)) / (p2(0) - p1(0));
    }
    if (std::abs(p2(1) - p1(1)) > eps) {
        a0y = (fGrid.y_min - p1(1)) / (p2(1) - p1(1));
        any = (fGrid.y_min + fGrid.dy * static_cast<MTfloat>(fGrid.ny) - p1(1)) / (p2(1) - p1(1));
    }
    if (std::abs(p2(2) - p1(2)) > eps) {
        a0z = (fGrid.z_min - p1(2)) / (p2(2) - p1(2));
        anz = (fGrid.z_min + fGrid.dz * static_cast<MTfloat>(fGrid.nz) - p1(2)) / (p2(2) - p1(2));
    }

    // calculate the minimum and maximum values of alpha 
    auto ammx = std::minmax(a0x, anx);
    auto ammy = std::minmax(a0y, any);
    auto ammz = std::minmax(a0z, anz);

    MTfloat amin = std::max(std::max(ammx.first, ammy.first), ammz.first);
    MTfloat amax = std::min(std::min(ammx.second, ammy.second), ammz.second);

    // return empty vector if the ray does not intersect with the grid
    if (amax <= amin) {
        return vpath;
    }

    // first intersect point alpha of x, y, z
    MTfloat ax = p1(0) + amin * (p2(0) - p1(0)))


}