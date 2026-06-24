#ifndef GRIDFLUIDSOLVER_H
#define GRIDFLUIDSOLVER_H

#include "GridManager.h"
#include "CenterScalarGrid.h"
#include "FluidSolver.h"
#include "SemiLagrangianAdvectSolver.h"
#include "PressureSolver.h"

class GridFluidSolver: public FluidSolver
{
public:
    GridFluidSolver(Size3 res, Vec3 origin = Vec3(0, 0, 0), Vec3 spacing = Vec3(1, 1, 1));
    virtual ~GridFluidSolver()
    {}

    Size3 res() const;
    Vec3 spacing() const;
    Vec3 origin() const;

public:
    GridManager &grids();
    void reinitializeFreeSurface();
    SemiLagrangianAdvectSolver advectionSolver;

    void resetSourceVelocity(double cylinderRadius,
        double cylinderHeight, 
        double offset,
        double sourceVelocity);

    std::shared_ptr<CenterScalarGrid> colliderSDF;

protected:
    void onAct(double time_step) override;

    virtual void calcExternalForces(double time_step);
    virtual void calcPressure(double time_step);
    virtual void calcAdvection(double time_step);

private:
    Size3 solverRes;
    Vec3 solverOrigin;
    Vec3 solverSpacing;
    Vec3 gravity = Vec3(0, -9.81, 0);
    GridManager solverGrids;
    //std::shared_ptr<CenterScalarGrid> colliderSDF;
    int count = 0;

    PressureSolver pressureSolver;

    void applyGravity(double time_step);
    void applyBoundaryCondition();

    void applyDomainBoundaryConditions(std::shared_ptr<MACVectorGrid> velocity);
    void applyColliderBoundaryConditions(std::shared_ptr<MACVectorGrid> velocity);

    void preAct(double time_step);
    void postAct(double time_step);

    void extrapolateIntoCollider(ScalarGrid &grid);
    void reinitializeFreeSurface(const ScalarGrid &inputMesh, double maxDistance, ScalarGrid &outputMesh);
    void getDerivatives(ScalarGrid &grid,
                        const Size3 indices,
                        std::array<double, 2> &dx,
                        std::array<double, 2> &dy,
                        std::array<double, 2> &dz);

    void extrapolateVelocityToAir();
};

#endif //GRIDFLUIDSOLVER_H
