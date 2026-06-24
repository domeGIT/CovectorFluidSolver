#ifndef PRESSURESOLVER_H
#define PRESSURESOLVER_H
#include "MACVectorGrid.h"
#include "ScalarGrid.h"

#include <eigen3/Eigen/Sparse>

enum Medium
{
    FLUID = 0,
    AIR = 1,
    COLLIDER = 2
};

class PressureSolver
{
public:
    PressureSolver();
    virtual ~PressureSolver();
    void solve(const MACVectorGrid &input,
               const ScalarGrid &colliderSDF,
               const ScalarGrid &fluidSDF,
               MACVectorGrid &output, double dt,
               const ScalarGrid* sourceMask = nullptr);

private:
    std::vector<Medium> markers;
    std::vector<char> sourceRegion;

    void setupMarkers(const MACVectorGrid &input, const ScalarGrid &collider, const ScalarGrid &fluid, const ScalarGrid* sourceMask = nullptr);
    void setupSystem(const MACVectorGrid &input, double dt);
    void applyPressureGradient(const MACVectorGrid &input, MACVectorGrid &output, double dt);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd x, b;
};
#endif //PRESSURESOLVER_H
