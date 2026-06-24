#ifndef VERTEXVECTORGRID_H
#define VERTEXVECTORGRID_H

#include "CollocatedGrid.h"
class VertexVectorGrid: public CollocatedGrid
{
public:
    Size3 dataSize() const override
    {
        return res() + Size3(1,1,1);
    }
    Eigen::Vector3d dataOrigin() const override
    {
        return origin();
    }
};


#endif //VERTEXVECTORGRID_H
