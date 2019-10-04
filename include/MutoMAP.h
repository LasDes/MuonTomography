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
    MTfloat fBeta;
    MTfloat fP;

    // some flags to enable addition features
    bool fUseAllNeighbours;
    bool fAddObjectPrior;


};