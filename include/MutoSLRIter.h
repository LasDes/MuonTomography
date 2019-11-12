/* Iterative version of SLR */

#pragma once

#include "MutoTypes.h"
#include "MutoSLR.h"
#include "MutoConfig.h"

class MutoSLRIter : public MutoSLR {

public:
    // default constructor, set default values for PoCA process
    MutoSLRIter(); 
    // construct with json configuration
    MutoSLRIter(json); 

    virtual ~MutoSLRIter();

    // implementation of reconstruction function usin SLR method
    virtual Image reconstruct(const MutoMuonData &);

private:
    int fNIteration;
};