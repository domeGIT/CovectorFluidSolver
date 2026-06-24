#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <eigen3/Eigen/Dense>
#include <array>
#include <vector>

#include <eigen3/Eigen/Dense>
typedef Eigen::Array<int, 3, 1> Size3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;

class BaseGrid;
class ScalarGrid;
class MACVectorGrid;

namespace Utils
{
/**
 * retieves the barycentic weight aint the ranges and a given value. return a value representing the an index
 * aint the range, and the weight given to this index. values above or below this ranged are clamped.
 * @param x the value aint this range. e.g. in [3,6) (keeping in mind that the upper range
 * is like an array iteration), 4.1 represents 1.1 units aint this range
 * @param low the lower boundary index
 * @param high the higher boundary index
 * @param index the caluclated index
 * @param weight the calculated weight
 */
inline void getBarycentric(const double x, int low, int high, int &index, double &weight)
{
    //floor the x given and pass as the index
    double s = std::floor(x);
    index = static_cast<int>(s);

    int oldLow = low;

    low = 0;
    high -= oldLow;

    if (low == high)
    {
        index = low;
        weight = 0;
    }

    else if (index < low)
    {
        index = low;
        weight = 0;
    }

    else if (index > high - 1)
    {
        index = high - 1;
        weight = 1;
    }
    else
    {
        weight = static_cast<double>(x - s);
    }

    index += oldLow;
}

void getCoordinatesAndWeights(const Size3 res,
                              const Vec3 spacing,
                              const Vec3 origin,
                              const Eigen::Vector3d x,
                              std::array<Size3, 8> &indices,
                              std::array<double, 8> &weights,
                              Vec3 *triweights = nullptr);

double trilinearInterpolateFromGrid(const ScalarGrid &grid, const Vec3 point);
Vec3 trilinearInterpolateFromGrid(const MACVectorGrid &grid, const Vec3 point);
double cubicInterpolateFromGrid(const ScalarGrid& grid, const Vec3 point);
Vec3 cubicInterpolateFromGrid(const MACVectorGrid& grid, const Vec3 point);

double catmullRom(std::array<double,4> values, double x);
double trilinearInterpolate(std::array<double, 8> values, Vec3 weights);
double bilinearInterpolate(std::array<double, 4> values, Vec2 weights);
double lerp(double x, double y, double weight);

void extrapolateToRegion(
    const std::vector<double> &input,
    const Size3 dims,
    const std::vector<bool> &validMatrix,
    unsigned int numIters,
    std::vector<double> &output);

double pseudoTimeStep(
    const ScalarGrid &sdf,
    const Vec3 gridSpacing,
    double maxCFL);
}

#endif //MATHUTILS_H
