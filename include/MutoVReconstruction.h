/* abstract class for future image reconstruction algorithms
   two virtual functions should be provided if extending this class
    + reconstruct : MutoMuonData -> Image
    + initialize : MutoConfig -> bool
*/

#pragma once

#include "MutoTypes.h"
#include "MutoConfig.h"

class MutoVReconstruction {
public:
    MutoVReconstruction() {}
    virtual ~MutoVReconstruction() = 0;

    virtual Image reconstruct(const MutoMuonData &) = 0;
    virtual bool initialize(const MutoConfig &) = 0;
};