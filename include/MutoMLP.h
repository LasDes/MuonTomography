/* Calculate most likely path of rays through voxels using Cubic fit method
    + initialize the class with voxel informations;
    + a function to return path through each voxel btw two points with directions (points should be at the top/bottom boundaries)
    + steps of MLP acquisition :
        - transform in/out rays and voxel such that incoming rays become "vertical";
        - fit a cubic line using 4 equations from in/out ray;
        - get interacting coordination with voxels in XZ and YZ plane respectively;
        - transform interacting points back into original voxel index.

   references: 
    Schulte, R. W., Penfold, S. N., Tafas, J. T., & Schubert, K. E. (2008). A maximum likelihood proton path formalism for application in proton computed tomography. Medical Physics, 35(11), 4849–4856. https://doi.org/10.1118/1.2986139
    Collins-Fekete, C. A., Volz, L., Portillo, S. K. N., Beaulieu, L., & Seco, J. (2017). A theoretical framework to predict the most likely ion path in particle imaging. Physics in Medicine and Biology, 62(5), 1777–1790. https://doi.org/10.1088/1361-6560/aa58ce
    Yi, H., Zeng, Z., Yu, B., Cheng, J., Zhao, Z., Wang, X., … Wang, Y. (2016). Bayesian-theory-based most probable trajectory reconstruction algorithm in cosmic ray muon tomography. 2014 IEEE Nuclear Science Symposium and Medical Imaging Conference, NSS/MIC 2014, 1–4. https://doi.org/10.1109/NSSMIC.2014.7431084
    Weihe's note: Explicit formula for Most Likely Path calculation
*/

#pragma once

#include "MutoTypes.h"
#include <vector>

class MutoMLP {

public:
    MutoMLP(const VoxelGrid &);
    ~MutoMLP() {}

    // get voxels btw in/out rays using simplified MLP algorithm
    // return values: vectors representing index and path length of each voxel
    std::vector<VoxelData> getVoxelPath(const Vector3 &, const Vector3 &);

private:
    VoxelGrid fGrid;
    std::vector<VoxelData> fPath;
};
