#ifndef VECTORGRIDBUILDER_H
#define VECTORGRIDBUILDER_H

#include "VectorGrid.h"
#include <memory>

class VectorGridBuilder
{
public:
    VectorGridBuilder()
    {}
    virtual ~VectorGridBuilder()
    {}
    virtual std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) = 0;

};

class VertexVectorBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};

class CenterVectorBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};

class MACGridBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};
#endif //VECTORGRIDBUILDER_H
