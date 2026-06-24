#ifndef MYVECTORFIELD_H
#define MYVECTORFIELD_H

#include "BaseGrid.h"

#include <eigen3/Eigen/Dense>
typedef Eigen::Array<int, 3, 1> Size3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;


class VectorGrid: public BaseGrid
{
public:
    VectorGrid()
    {};
    virtual ~VectorGrid()
    {};


    virtual void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing) override
    {
        BaseGrid::setSize(res, origin, spacing);
        onResize(res, spacing, origin, Vec3(0,0,0));
    }

    void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing, Vec3 initValue)
    {
        BaseGrid::setSize(res, origin, spacing);
        onResize(res, spacing, origin, initValue);
    }

    virtual Vec3 sample(const Vec3 point) = 0;

protected:
    virtual void onResize(const Size3& resolution, Eigen::Vector3d gridSpacing,
                          Eigen::Vector3d origin, Eigen::Vector3d inittialValue) = 0;


};
#endif //MYVECTORFIELD_H
