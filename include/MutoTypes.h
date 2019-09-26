/* This file defines all useful types for the entire program. */

#pragma once

// define floating point in the whole program (reserved for GPU)
typedef double MTfloat;
typedef unsigned int MTindex;

// use Matrix library: eigen
#include <Eigen/Dense>
#include <Eigen_unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;

typedef Matrix<MTfloat, 2, 2> Mat2x2;
typedef Matrix<MTfloat, 1, 3> Vector3;

// image class to store reconstructed scattering density map
typedef Tensor<MTfloat, 3> Image;

// voxel grid 
struct VoxelGrid {
	MTindex nx;
	MTindex ny;
	MTindex nz;

	MTfloat x_min;
	MTfloat y_min;
	MTfloat z_min;

	MTfloat dx;
	MTfloat dy;
	MTfloat dz;
};

// voxel index and path length
struct VoxelData {
    MTindex x;
    MTindex y;
    MTindex z;
    MTfloat length;
};