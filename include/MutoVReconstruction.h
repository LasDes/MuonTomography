/* abstract class for future image reconstruction algorithms
   two virtual functions should be provided if extending this class
    + reconstruct : MutoMuonData -> Image
    + initialize : MutoConfig -> bool
*/

#pragma once

#include "MutoTypes.h"

class MutoVReconstruction {
public:
    MutoVReconstruction() {}
    virtual ~MutoVReconstruction() {}

    virtual Image reconstruct(const MutoMuonData &) = 0;
};