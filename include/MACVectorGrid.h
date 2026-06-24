#ifndef MACVECTORGRID_H
#define MACVECTORGRID_H
#include <iostream>
#include "CollocatedGrid.h"
class MACVectorGrid: public VectorGrid
{
public :
    MACVectorGrid() {};

    MACVectorGrid(const MACVectorGrid &other)
    {
        setSize(other.res(), other.origin(), other.spacing());

        for (unsigned int i = 0; i < other.dataU.size(); i++)
            this->dataU[i] = other.dataU[i];

        for (unsigned int i = 0; i < other.dataV.size(); i++)
            this->dataV[i] = other.dataV[i];

        for (unsigned int i = 0; i < other.dataW.size(); i++)
            this->dataW[i] = other.dataW[i];
    }


    const double u(int i, int j, int k) const;
    const double v(int i, int j, int k) const;
    const double w(int i, int j, int k) const;

    double u(int i, int j, int k);
    double v(int i, int j, int k);
    double w(int i, int j, int k);

    double divergenceAt(int i, int j, int k) const;
    double divergenceAt(Size3 indices) const;
    double divergence(Vec3 x) const;
    double divergence(int i, int j, int k) const;

    Vec3 curlAt(Size3 indices) const;
    Vec3 curlAt(int i, int j, int k) const;
    Vec3 curl(const Vec3 x) const;

    double maxNorm() const;

    void fillData(std::function<Vec3(Vec3)> predicate);
    void fillData(std::function<Vec3(double, double, double)> predicate);
    void fillData(double value);

    void setUDataAt(int i, int j, int k, double value);
    void setVDataAt(int i, int j, int k, double value);
    void setWDataAt(int i, int j, int k, double value);

    void setUDataAt(Size3 indices, double value);
    void setVDataAt(Size3 indices, double value);
    void setWDataAt(Size3 indices, double value);

    void getDataU(std::vector<double>& result) const;
    void getDataV(std::vector<double>& result) const;
    void getDataW(std::vector<double>& result) const;


    Vec3 dataCenterPosU(int i, int j, int k) const;
    Vec3 dataCenterPosU(const Vec3 indices) const;

    Vec3 dataCenterPosV(int i, int j, int k) const;
    Vec3 dataCenterPosV(const Vec3 indices) const;

    Vec3 dataCenterPosW(int i, int j, int k) const;
    Vec3 dataCenterPosW(const Vec3 indices) const;

    void forEachU(std::function<void(int, int, int)> operation);
    void forEachV(std::function<void(int, int, int)> operation);
    void forEachW(std::function<void(int, int, int)> operation);

    double sizeU();
    double sizeV();
    double sizeW();

    Vec3 atCenter(int i, int j, int k) const;

    Vec3 sample(const Vec3 point) override;
    Vec3 cubicSample(const Vec3 point);

    inline int flatIndexU(int i, int j, int k) const
    { return i * getUDims()[1] * getUDims()[2] + j * getUDims()[2] + k; }
    inline int flatIndexV(int i, int j, int k) const
    { return i * getVDims()[1] * getVDims()[2] + j * getVDims()[2] + k; }
    inline int flatIndexW(int i, int j, int k) const
    { return i * getWDims()[1] * getWDims()[2] + j * getWDims()[2] + k;}

    inline Size3 getUDims() const
    { return Size3(res()[0] + 1, res()[1], res()[2]); }
    inline Size3 getVDims() const
    { return Size3(res()[0], res()[1] + 1, res()[2]); }
    inline Size3 getWDims() const
    { return Size3(res()[0], res()[1], res()[2] + 1); }

    inline Vec3 getUOrigin() const { return originU; }
    inline Vec3 getVOrigin() const { return originV; }
    inline Vec3 getWOrigin() const { return originW; }

#if 0
    ///       (i, j, k)
    double vorticityMagnitude(int i, int j, int k) const;

    ///       (i, j, k)
    Vec3 vorticityVector(int i, int j, int k) const;

    ///      
    double computeKineticEnergy() const;

    ///     
    double computeEnstrophy() const;

    ///     ( CFL  .)
    double maxVorticity() const;
#endif

// Interpolacija pojedinačnih komponenti
double sampleUAt(const Vec3 point) const;
double sampleVAt(const Vec3 point) const;
double sampleWAt(const Vec3 point) const;

// Interpolacija celokupne brzine
Vec3 sampleVelocityAt(const Vec3 point) const;

protected:

    void onResize(const Size3 &resolution,
                  Vec3 gridSpacing,
                  Vec3 origin,
                  Eigen::Vector3d initialValue) override;
private:

    std::vector<double> dataU;
    std::vector<double> dataV;
    std::vector<double> dataW;

    Vec3 originU,originV, originW;
};
#endif //MACVECTORGRID_H
