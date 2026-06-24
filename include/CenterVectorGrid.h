#ifndef CENTERVECTORGRID_H
#define CENTERVECTORGRID_H

#include "CollocatedGrid.h"
class CenterVectorGrid: public CollocatedGrid
{
public:
    Size3 dataSize() const override
    {
        return res();
    }
    Eigen::Vector3d dataOrigin() const override
    {
        return origin() + 0.5 * spacing();
    }
};

#endif //CENTERVECTORGRID_H
