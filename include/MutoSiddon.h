/* Calculate ray (straight) intersection with voxel using modified Siddon Algorithm
    + initialize the class with voxel informations;
    + a function to return path through each voxel btw two points

   reference: 
    Jacobs, F. (1998). A fast algorithm to calculate the exact radiological path through a pixel or voxel space.
*/

#pragma once

#include "MutoTypes.h"
#include <vector>

class MutoSiddon {

public:
    MutoSiddon(const VoxelGrid &);
    ~MutoSiddon() {;}

    // get voxels btw two points using Siddon algorithm
    // return values: vectors representing index and path length of each voxel
    std::vector<VoxelData> getVoxelPath(const Vector3 &, const Vector3 &);

private:
    VoxelGrid fGrid;
};
