/* add penalty to ML and optimize using surrogate function, a.k.a Maximum a posterior method
    + load initial image e.g. PoCA, or random generated noise;
    + optimize ML part using the same method in MLSD/EM algorithm;
    + penalize the object function with l-norm of neighbouring voxels;
    + optimize the penalty part by solving a 3rd degree polynominal equation.

   reference: 
    Yu, B., Zhao, Z., Wang, X., Wu, D., Zeng, Z., Zeng, M., … Cheng, J. (2016). A unified framework for penalized statistical muon tomography reconstruction with edge preservation priors of lp norm type. Nuclear Instruments and Methods in Physics Research, Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 806, 199–205. https://doi.org/10.1016/j.nima.2015.09.113
    Zhang, H., Wang, J., Zeng, D., Tao, X., & Ma, J. (2018). Regularization strategies in statistical image reconstruction of low-dose x-ray CT: A review. Medical Physics, 45(10), e886–e907. https://doi.org/10.1002/mp.13123
*/

#pragma once

#include "MutoMLSD.h"

class MutoMAP : public MutoMLSD {

public:
    // default constructor, set default values for MLSD process
    MutoMAP(); 
    // construct with json configuration
    MutoMAP(json); 

    virtual ~MutoMAP();

    // implementation of reconstruction function usin MAP method
    virtual Image reconstruct(const MutoMuonData &);

private:
    // private members specific for regularization in MAP algorithm
    //  + beta: strength of regularization penalty
    //  + p: lp-norm of the pair-wise potential function
    //  + delta: a small number that keep function derivable at 0
    MTfloat fBeta;
    MTfloat fP;
    MTfloat fDelta;

    // some flags to enable addition features
    //  + all_neighbours: instead of 6 direct neighbours, use all 26 blocks (weighted by dist.) around the center voxel
    //  + add_object_prior: add independent prior of center voxel
    bool fUseAllNeighbours;
    bool fAddObjectPrior;

    // not configurable internal parameters
    MTfloat fUnitConv; // convert from rad^2/mm to mrad^2/cm

    // default initialization
    void defaultInitialize();

protected:
    // calculate coefficients from neighbouring voxels
    std::pair<std::vector<MTfloat>, std::vector<MTfloat>> calculateCoefFromNeighbours(MTindex j, const Image& img, MTindex Mj, 
                                                                                      MTfloat& p2, MTfloat& p3);
    Matrix<MTindex, 3, 1> index2xyz(MTindex i);
    std::vector<MTfloat> getSixNeighbours(MTindex j, const Image& img);
    std::vector<MTfloat> getAllNeighbours(MTindex j, const Image& img);
    MTfloat chooseBestEstimation(MTindex Mj, MTfloat lj, MTfloat lj_ML, MTfloat* roots, MTindex nRoots, 
                                 const std::pair<std::vector<MTfloat>, std::vector<MTfloat>>& weightAndLambda);

    MTfloat surrogateFunctionValue(MTindex Mj, MTfloat ljn, MTfloat ljML, MTfloat lj, std::vector<MTfloat> w, std::vector<MTfloat> lm);

    // cubic equation solver using Newton method (x10 faster than PolynomialSolver in Eigen)
    MTindex cubicEqnSolver(MTfloat* roots, MTfloat p3, MTfloat p2, MTfloat p1, MTfloat p0, MTfloat tol);
};