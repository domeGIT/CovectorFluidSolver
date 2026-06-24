#ifndef SCALARGRIDBUILDER_H
#define SCALARGRIDBUILDER_H

#include <ScalarGrid.h>
#include <memory>

class ScalarGridBuilder
{
public:
    ScalarGridBuilder()
    {}
    virtual ~ScalarGridBuilder()
    {}
    virtual std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) = 0;

};

class VertexScalarBuilder: public ScalarGridBuilder
{
public:
    std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) override;

};

class CenterScalarBuilder: public ScalarGridBuilder
{
public:

    std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) override;

};
#endif //SCALARGRIDBUILDER_H
