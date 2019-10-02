#include "MutoSiddon.h"
#include <cmath>

MutoSiddon::MutoSiddon(const VoxelGrid & grid) {
    fGrid = grid;
    fPath.reserve(2*fGrid.nz);
}

std::vector<VoxelData> MutoSiddon::getVoxelPath(const Vector3 & p1, const Vector3 & p2) {
    // initialize return data
    fPath.clear();

    // Enclidean distance btw two points
    MTfloat length = std::sqrt((p2 - p1).dot(p2 - p1));

    // get first and last intersection points for each parallel planes
    // i.e. amin and amax for x, y, z 
    MTfloat a0x = -0.1 , a0y = -0.1, a0z = -0.1, anx = 1.1, any = 1.1, anz = 1.1;
    MTfloat eps = 1.0e-6;
    MTfloat dp21x = p2(0) - p1(0);
    MTfloat dp21y = p2(1) - p1(1);
    MTfloat dp21z = p2(2) - p1(2);


    if (eps < std::abs(dp21x)) { // check for parallel rays
        a0x = (fGrid.x_min - p1(0)) / dp21x;
        anx = (fGrid.x_min + fGrid.dx * static_cast<MTfloat>(fGrid.nx) - p1(0)) / dp21x;
    } else if (p1(0) < fGrid.x_min || fGrid.x_min + fGrid.dx*static_cast<MTfloat>(fGrid.nx) <= p1(0)) { 
        // return empty vector if parallel ray is out of voxel world
        return fPath;
    }
    if (eps < std::abs(dp21y)) {
        a0y = (fGrid.y_min - p1(1)) / dp21y;
        any = (fGrid.y_min + fGrid.dy * static_cast<MTfloat>(fGrid.ny) - p1(1)) / dp21y;
    } else if (p1(1) < fGrid.y_min || fGrid.y_min + fGrid.dy*static_cast<MTfloat>(fGrid.ny) <= p1(1)) {
        return fPath;
    }
    if (eps < std::abs(dp21z)) {
        a0z = (fGrid.z_min - p1(2)) / dp21z;
        anz = (fGrid.z_min + fGrid.dz * static_cast<MTfloat>(fGrid.nz) - p1(2)) / dp21z;
    } else if (p1(2) < fGrid.z_min || fGrid.z_min + fGrid.dz*static_cast<MTfloat>(fGrid.nz) <= p1(2)) {
        return fPath;
    }

    // calculate the minimum and maximum values of alpha 
    auto ammx = std::minmax(a0x, anx);
    auto ammy = std::minmax(a0y, any);
    auto ammz = std::minmax(a0z, anz);

    MTfloat amin = std::max(std::max(std::max(ammx.first, ammy.first), ammz.first), 0.0);
    MTfloat amax = std::min(std::min(std::min(ammx.second, ammy.second), ammz.second), 1.0);

    // return empty vector if the ray does not intersect with the grid
    if (amax <= amin) {
        return fPath;
    }

    // first intersect point alpha of x, y, z
    // initialize current alpha to its minimum value
    MTfloat alpha = amin; 

    MTindex ix = static_cast<MTindex>((alpha*dp21x + p1(0) - fGrid.x_min) / fGrid.dx + (dp21x <= 0 ? -eps : eps)); // add/minus eps 
    MTindex iy = static_cast<MTindex>((alpha*dp21y + p1(1) - fGrid.y_min) / fGrid.dy + (dp21y <= 0 ? -eps : eps));
    MTindex iz = static_cast<MTindex>((alpha*dp21z + p1(2) - fGrid.z_min) / fGrid.dz + (dp21z <= 0 ? -eps : eps));

    // updating parameters for x, y, z in in the loop
    int dix = 0, diy = 0, diz = 0;    // initialize to 0 for rays that are parallel 
    MTfloat dax = 0.0, day = 0.0, daz = 0.0;

    if (eps < std::abs(dp21x)) {
        dax = fGrid.dx / (std::abs(dp21x));
        dix = dp21x < 0 ? -1 : 1;
    }
    if (eps < std::abs(dp21y)) {
        day = fGrid.dy / (std::abs(dp21y));
        diy = dp21y < 0 ? -1 : 1;
    }
    if (eps < std::abs(dp21z)) {
        daz = fGrid.dz / (std::abs(dp21z));
        diz = dp21z < 0 ? -1 : 1;
    }

    // ax, ay, az for next intersect plane after the first voxel
    // (!this part of code is quite tricky and not easy to understand)
    MTfloat ax = 1.1, ay = 1.1, az = 1.1; // initialize to some value larger than 1.0

    if (std::abs(amin - ammx.first) < eps) { // first point == alpha_min 
        // (NOTE: alphas for parallel ray were initialized out of the range)
        ax = amin + dax;
    } else if (dix != 0) { // ray not parallel with plane
        // add dax if index is increasing
        ax = (static_cast<MTfloat>(ix) * fGrid.dx + fGrid.x_min - p1(0)) / dp21x + (dix == 1 ? dax : 0.0);
    }
    if (std::abs(amin - ammy.first) < eps) {
        ay = amin + day;
    } else if (diy != 0)  {
        ay = (static_cast<MTfloat>(iy) * fGrid.dy + fGrid.y_min - p1(1)) / dp21y + (diy == 1 ? day : 0.0);
    }
    if (std::abs(amin - ammz.first) < eps) {
        az = amin + daz;
    } else if (diz != 0) {
        az = (static_cast<MTfloat>(iz) * fGrid.dz + fGrid.z_min - p1(2)) / dp21z + (diz == 1 ? daz : 0.0);
    }

    // store previous alpha value
    MTfloat alpha_last = alpha;
    // loop through the ray and the voxel grid
    while (alpha < amax) {
        // current voxel index
        VoxelData voxel;
        voxel.x = ix;
        voxel.y = iy;
        voxel.z = iz;

        // get minimum alpha (i.e. current intersect plane) among ax, ay, az
        if (ax < ay && ax < az) { // ax is current minimum 
            ix += dix;
            alpha = ax;
            ax += dax;
        } else if (ay < az) { // ay is current minimum
            iy += diy;
            alpha = ay;
            ay += day;
        } else { // az is current minimum
            iz += diz;
            alpha = az;
            az += daz;
        }

        // get length of ray in each voxel and push data to output vector
        MTfloat dalpha = std::min(alpha, amax) - alpha_last;
        if (eps < dalpha) {
            voxel.length = dalpha * length;
            fPath.push_back(voxel);
            // refreshs alpha_last value if length is recorded (corner values are added to next voxel)
            alpha_last = alpha;        
        }
    }
    
    return fPath;
}