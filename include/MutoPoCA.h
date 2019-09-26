/* Use Point of Closest Approach (PoCA) method to get muon image
    + get the closest point btw the in/out ray as scattering point;
    + distribute weight to each passing voxel;
    + summing weight for all the recorded muons; 

   reference: 
    Schultz, L. J. O. E. (2003). COSMIC RAY MUON RADIOGRAPHY. PhD thesis
*/

#pragma once

#include "MutoTypes.h"
#include "MutoVReconstruction.h"

class MutoPoCA : public MutoVReconstruction {

};