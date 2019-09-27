/* This file defines all useful types for the entire program. */

#pragma once

#include <vector>
#include <string>

// define floating point in the whole program (reserved for GPU)
typedef double MTfloat;
typedef unsigned int MTindex;

// use Matrix library: eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

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

// muon ray data containing energy (future), position and direction
struct RayData {
	MTfloat ein;
	MTfloat eout;
	Vector3 pin;
	Vector3 pout;
	Vector3 din;
	Vector3 dout;
};

typedef std::vector<RayData> MutoMuonData;

struct ImageHeader {
	std::string method; // only the first 16 chars will be write to file header  
	VoxelGrid grid;
	MTindex nRay;
	MTindex nIter;
	MTindex fPrecision;
};