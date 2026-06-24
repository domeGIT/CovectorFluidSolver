#ifndef BaseGrid_H
#define BaseGrid_H

#include <eigen3/Eigen/Dense>

#include <eigen3/Eigen/Dense>
typedef Eigen::Array<int, 3, 1> Size3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;


class BoundingBox
{
    Eigen::Vector3d min, max;

public:
    BoundingBox(const Eigen::Vector3d min, const Eigen::Vector3d max)
        : min(min), max(max)
    {

    }

    BoundingBox()
        : min(Vec3(0, 0, 0)), max(Vec3(0, 0, 0))
    {

    }

    const Eigen::Vector3d& getMin() const
    {
        return min;
    }

    const Eigen::Vector3d& getMax() const
    {
        return max;
    }
};

class BaseGrid
{
public:

    BaseGrid(){};
    virtual ~BaseGrid(){};

    inline const Size3 res() const{ return gridRes; };
    inline const Vec3 origin() const{ return gridOrigin;};
    inline const Vec3 spacing() const{ return gridSpacing;};
    inline const BoundingBox getBounds() const
    { return bounds; };

    Vec3 cellCenterPosition(Size3 indices) const
    {
        Vec3 centerIndices(indices[0] + 0.5,indices[1] + 0.5,indices[2] + 0.5);
        Vec3 offsets = centerIndices.cwiseProduct(gridSpacing);
        return gridOrigin + offsets;
    };

    virtual void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing)
    {
        this->gridRes = res;
        this->gridOrigin = origin;
        this->gridSpacing = spacing;

        Vec3 dGridRes(gridRes[0],gridRes[1],gridRes[2]);
        dGridRes = dGridRes.cwiseProduct(gridSpacing);

        bounds = BoundingBox(origin, origin + dGridRes);
    }

    void setSizeDefault(Size3 res)
    {
        setSize(res,Vec3(0,0,0),Vec3(1,1,1));
    }

private:
    Size3 gridRes;
    Vec3 gridSpacing = Vec3(1, 1, 1);
    
    Vec3 gridOrigin = Vec3(0,0,0);
    BoundingBox bounds = BoundingBox(Vec3(1, 1, 1), Vec3(-1, -1, -1));

};

#endif //BaseGrid_H
