#include <memory>
#include <VertexVectorGrid.h>
#include <CenterVectorGrid.h>
#include <MACVectorGrid.h>
#include "VectorBuilders.h"

#include <eigen3/Eigen/Dense>
typedef Eigen::Array<int, 3, 1> Size3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;


std::shared_ptr<VectorGrid> VertexVectorBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<VertexVectorGrid> grid = std::make_shared<VertexVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}

std::shared_ptr<VectorGrid> CenterVectorBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<CenterVectorGrid> grid = std::make_shared<CenterVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}

std::shared_ptr<VectorGrid> MACGridBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<MACVectorGrid> grid = std::make_shared<MACVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}
