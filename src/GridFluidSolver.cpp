#include <Utils.h>
#include <iomanip>
#include "GridFluidSolver.h"

extern bool USE_COVECTOR_CORRECTION;
extern bool USE_COLLIDER;
extern int Scenario;
extern double cylinderRadius;
extern double cylinderHeight;
extern double cylinderOffset;
extern double sourceVelocity;

double subDt;

GridFluidSolver::GridFluidSolver(Size3 res, Vec3 origin, Vec3 spacing)
    :
    solverRes(res), solverOrigin(origin), solverSpacing(spacing), solverGrids(res, origin, spacing)
{
    CenterScalarBuilder free_surfaceBuilder;
    solverGrids.addScalarGrid("FREE_SURFACE", free_surfaceBuilder);
    solverGrids.addScalarGrid("VORTICITY", free_surfaceBuilder);
    solverGrids.addScalarGrid("SOURCE_MASK", free_surfaceBuilder);

    colliderSDF = std::make_shared<CenterScalarGrid>(res, origin, spacing, std::numeric_limits<double>().max());

    if(USE_COLLIDER) // ako postoji kolajder potrebno ga je ovde definisati!

    colliderSDF->fillData([&](double x, double y, double z) -> double
        {
           if(Scenario == 1)
           {
            // 1. Sferni kolajder: rastojanje od take do povrine sfere

            double d_x = x - 0.5 * spacing[0] * res[0];  // centar sfere (2, 1.12, 2)
            double d_y = y - 0.3 * spacing[1] * res[1];
            double d_z = z - 0.5 * spacing[2] * res[2];
            double distToCenter = std::sqrt(d_x * d_x + d_y * d_y + d_z * d_z);
            double sdfSphere = distToCenter - 0.0625 * spacing[1] * res[1]; // poluprecnik sfere 0.25
            return sdfSphere;
           }
        });
}

void GridFluidSolver::onAct(double time_step)
{
    preAct(time_step);

    calcExternalForces(time_step);
    calcAdvection(time_step);
    calcPressure(time_step);

    postAct(time_step);
}

void GridFluidSolver::calcExternalForces(double time_step)
{
    std::cout << "calcExternalForces ..." << std::endl;
    std::cout.flush();

    applyGravity(time_step);

    applyBoundaryCondition();

    if (Scenario == 6) {
        resetSourceVelocity(cylinderRadius, cylinderHeight, cylinderOffset, sourceVelocity);
    }

}

void GridFluidSolver::calcAdvection(double delta)
{
    std::cout << "calcAdvection ..." << std::endl;
    std::cout.flush();

    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("FREE_SURFACE"));
    auto sdfInput = std::make_shared<CenterScalarGrid>(*sdf);
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());

    extrapolateVelocityToAir();

    std::shared_ptr<MACVectorGrid> velocityPrev = std::make_shared<MACVectorGrid>(*velocityGrid);

    // ========================================================================
    // CFL-BASED SUB-STEPPING (Section 5.4.3)
    // ========================================================================

    double maxVelocity = 0.0;
    velocityPrev->forEachU([&](int i, int j, int k) {
        maxVelocity = std::max(maxVelocity, std::abs(velocityPrev->u(i, j, k)));
        });
    velocityPrev->forEachV([&](int i, int j, int k) {
        maxVelocity = std::max(maxVelocity, std::abs(velocityPrev->v(i, j, k)));
        });
    velocityPrev->forEachW([&](int i, int j, int k) {
        maxVelocity = std::max(maxVelocity, std::abs(velocityPrev->w(i, j, k)));
        });

    double dx = velocityGrid->spacing()[0];
    double cflLimit = 0.5;  // Preporuka iz rada (Section 5.4.3)
    int numSubSteps = 1;

    if (maxVelocity > 1e-8) {
        numSubSteps = std::max(1, (int)std::ceil(maxVelocity * delta / (dx * cflLimit)));
    }

    // Ogranii max substepova da ne uspori previe
    numSubSteps = std::min(numSubSteps, 10);

    subDt = delta / numSubSteps;

    std::cout << "CFL substepping: " << numSubSteps << " steps, dt=" << subDt
        << ", maxVel=" << maxVelocity << std::endl;
    std::cout.flush();

    // ========================================================================
    // Advekcija kroz sub-stepove
    // ========================================================================

    for (int step = 0; step < numSubSteps; ++step) {
        // 1. Advekcija LS funkcije
        advectionSolver.advect(*sdfInput, velocityPrev, subDt, *sdf, colliderSDF);

        // 2. Standardna SL advekcija brzine
        advectionSolver.advect(*velocityPrev, velocityPrev, subDt, *velocityGrid, colliderSDF);

        // 3. Kovektorska korekcija (ako je ukljuena)
        if (USE_COVECTOR_CORRECTION) {
            auto velocity_SL = std::make_shared<MACVectorGrid>(*velocityGrid);

            std::cout << "covector correction ..." << std::endl;

            advectionSolver.covector_advection_1step(
                velocityPrev, subDt, velocity_SL, velocityGrid, colliderSDF);
        }

        applyBoundaryCondition();

    // Re-enforce face-centered source velocity after BC in case BC overwrote faces
    if (Scenario == 6) {
        resetSourceVelocity(cylinderRadius, cylinderHeight, cylinderOffset, sourceVelocity);
    }

        // Auriraj velocityPrev za sledeci sub-step
        if (step < numSubSteps - 1) {
            velocityPrev = std::make_shared<MACVectorGrid>(*velocityGrid);
        }
    }
}

void GridFluidSolver::calcPressure(double time_step)
{
    std::cout << "calcPressure ..." << std::endl;
    std::cout.flush();

    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto velInput = std::make_shared<MACVectorGrid>(*velocity);
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("FREE_SURFACE"));

    // For Scenario 6 create a source mask marking the two cylindrical sources so pressure solver enforces Dirichlet inflow

    CenterScalarGrid *sourceMaskPtr = nullptr;
    std::unique_ptr<CenterScalarGrid> sourceMaskOwner;
    if (Scenario == 6)
    {
        auto sourceMask = std::make_unique<CenterScalarGrid>(*sdf);
        // initialize to positive (no source)
        sourceMask->fillData([](double x, double y, double z) { return 5.0; });

        Size3 res = sourceMask->dataSize();
        Vec3 spacing = sourceMask->spacing();

        double Lx = res[0] * spacing[0];
        double Ly = res[1] * spacing[1];
        double Lz = res[2] * spacing[2];

        Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);
        Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);

        for (int i = 0; i < res[0]; ++i) {
            for (int j = 0; j < res[1]; ++j) {
                for (int k = 0; k < res[2]; ++k) {
                    double x = (i + 0.5) * spacing[0];
                    double y = (j + 0.5) * spacing[1];
                    double z = (k + 0.5) * spacing[2];

                    double dy1 = y - cyl1_center[1];
                    double dz1 = z - cyl1_center[2];
                    double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
                    double axial_dist1 = x - cylinderOffset;
                    bool inside1 = (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight);

                    double dx2 = x - cyl2_center[0];
                    double dy2 = y - cyl2_center[1];
                    double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
                    double axial_dist2 = (Lz - cylinderOffset) - z;
                    bool inside2 = (radial_dist2 < cylinderRadius && axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight);

                    if (inside1 || inside2)
                        sourceMask->setDataAt(Size3(i, j, k), -1.0);
                    else
                        sourceMask->setDataAt(Size3(i, j, k), 5.0);
                }
            }
        }

        // Also mark the solver's FREE_SURFACE SDF as fluid inside the sources so extrapolation/markers treat them as fluid
        for (int i = 0; i < res[0]; ++i) {
            for (int j = 0; j < res[1]; ++j) {
                for (int k = 0; k < res[2]; ++k) {
                    if (sourceMask->at(i,j,k) < 0.0)
                        sdf->setDataAt(Size3(i,j,k), -1.0);
                }
            }
        }

        sourceMaskPtr = sourceMask.get();
        sourceMaskOwner = std::move(sourceMask);
    }

    // Diagnostics: count source faces and report max velocities before solve

    if (Scenario == 6 && sourceMaskPtr != nullptr)
    {
        int uCount = 0, wCount = 0;
        double uMaxBefore = 0.0, wMaxBefore = 0.0;
        Size3 uDims = velocity->getUDims();
        Size3 wDims = velocity->getWDims();
        Vec3 spacingGrid = solverGrids.spacing();
        double Lx = solverGrids.res()[0] * spacingGrid[0];
        double Ly = solverGrids.res()[1] * spacingGrid[1];
        double Lz = solverGrids.res()[2] * spacingGrid[2];
        Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);
        Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);

        for (int i = 0; i < uDims[0]; ++i)
            for (int j = 0; j < uDims[1]; ++j)
                for (int k = 0; k < uDims[2]; ++k)
                {
                    Vec3 pos = velocity->dataCenterPosU(i, j, k);
                    double x = pos[0], y = pos[1], z = pos[2];
                    double dy1 = y - cyl1_center[1];
                    double dz1 = z - cyl1_center[2];
                    double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
                    double axial_dist1 = x - cylinderOffset;
                    if (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight)
                    {
                        ++uCount;
                        uMaxBefore = std::max(uMaxBefore, velInput->u(i,j,k));
                    }
                }

        for (int i = 0; i < wDims[0]; ++i)
            for (int j = 0; j < wDims[1]; ++j)
                for (int k = 0; k < wDims[2]; ++k)
                {
                    Vec3 pos = velocity->dataCenterPosW(i, j, k);
                    double x = pos[0], y = pos[1], z = pos[2];
                    double dx2 = x - cyl2_center[0];
                    double dy2 = y - cyl2_center[1];
                    double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
                    double axial_dist2 = (Lz - cylinderOffset) - z;
                    if (radial_dist2 < cylinderRadius && axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight)
                    {
                        ++wCount;
                        wMaxBefore = std::max(wMaxBefore, velInput->w(i,j,k));
                    }
                }

        //std::cout << "Scenario6: protected source faces U=" << uCount << " W=" << wCount
            //<< " maxBefore U=" << uMaxBefore << " W=" << wMaxBefore << std::endl;
    }

    pressureSolver.solve(*velInput, *colliderSDF, *sdf, *velocity, time_step, sourceMaskPtr);

    applyBoundaryCondition();

    // Diagnostics: report max velocities after solve in source faces
    if (Scenario == 6 && sourceMaskPtr != nullptr)
    {
        int uCount = 0, wCount = 0;
        double uMaxAfter = 0.0, wMaxAfter = 0.0;
        Size3 uDims = velocity->getUDims();
        Size3 wDims = velocity->getWDims();
        Vec3 spacingGrid = solverGrids.spacing();
        double Lx = solverGrids.res()[0] * spacingGrid[0];
        double Ly = solverGrids.res()[1] * spacingGrid[1];
        double Lz = solverGrids.res()[2] * spacingGrid[2];
        Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);
        Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);

        for (int i = 0; i < uDims[0]; ++i)
            for (int j = 0; j < uDims[1]; ++j)
                for (int k = 0; k < uDims[2]; ++k)
                {
                    Vec3 pos = velocity->dataCenterPosU(i, j, k);
                    double x = pos[0], y = pos[1], z = pos[2];
                    double dy1 = y - cyl1_center[1];
                    double dz1 = z - cyl1_center[2];
                    double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
                    double axial_dist1 = x - cylinderOffset;
                    if (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight)
                    {
                        ++uCount;
                        uMaxAfter = std::max(uMaxAfter, velocity->u(i,j,k));
                    }
                }

        for (int i = 0; i < wDims[0]; ++i)
            for (int j = 0; j < wDims[1]; ++j)
                for (int k = 0; k < wDims[2]; ++k)
                {
                    Vec3 pos = velocity->dataCenterPosW(i, j, k);
                    double x = pos[0], y = pos[1], z = pos[2];
                    double dx2 = x - cyl2_center[0];
                    double dy2 = y - cyl2_center[1];
                    double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
                    double axial_dist2 = (Lz - cylinderOffset) - z;
                    if (radial_dist2 < cylinderRadius && axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight)
                    {
                        ++wCount;
                        wMaxAfter = std::max(wMaxAfter, velocity->w(i,j,k));
                    }
                }

        //std::cout << "Scenario6: protected source faces after solve U=" << uCount << " W=" << wCount
            //<< " maxAfter U=" << uMaxAfter << " W=" << wMaxAfter << std::endl;
    }

    // Re-enforce persistent source velocities for Scenario 6 after projection and boundary conditions
    if (Scenario == 6)
    {
        resetSourceVelocity(cylinderRadius, cylinderHeight, cylinderOffset, sourceVelocity);
    }

}

void GridFluidSolver::preAct(double time_step)
{
    assert(grids().spacing()[0] == grids().getVelocityGrid()->spacing()[0]);
}

void GridFluidSolver::postAct(double time_step)
{
    reinitializeFreeSurface();
}
void GridFluidSolver::applyGravity(double time_step)
{
    std::cout << "applyGravity ..." << std::endl;

    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    assert(velocityGrid != nullptr);

    velocityGrid->forEachU([&](int i, int j, int k)
                           {
                               double uValue = velocityGrid->u(i, j, k);
                               velocityGrid->setUDataAt(i, j, k, uValue + time_step * gravity[0]);
                           });

    velocityGrid->forEachV([&](int i, int j, int k)
                           {
                               double vValue = velocityGrid->v(i, j, k);
                               velocityGrid->setVDataAt(i, j, k, vValue + time_step * gravity[1]);
                           });

    velocityGrid->forEachW([&](int i, int j, int k)
                           {
                               double wValue = velocityGrid->w(i, j, k);
                               velocityGrid->setWDataAt(i, j, k, wValue + time_step * gravity[2]);
                           });
}

void GridFluidSolver::applyBoundaryCondition()
{
    std::cout << "applyBoundaryCondition ..." << std::endl;

    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    assert(velocity != nullptr);

    // Prvo primenimo granine uslove na granicama domena
    applyDomainBoundaryConditions(velocity);

    // Zatim primenimo granine uslove na kolajderu ako postoji
    if (USE_COLLIDER) {
        applyColliderBoundaryConditions(velocity);
    }
}

void GridFluidSolver::applyDomainBoundaryConditions(std::shared_ptr<MACVectorGrid> velocity)
{
    Size3 uDims = velocity->getUDims();
    Size3 vDims = velocity->getVDims();
    Size3 wDims = velocity->getWDims();

    // Free-slip na granicama domena: normalna komponenta = 0

    // X pravac (U komponenta) - leva i desna granica
    for (int j = 0; j < uDims[1]; j++)
    {
        for (int k = 0; k < uDims[2]; k++)
        {
            velocity->setUDataAt(0, j, k, 0);          // leva granica
            velocity->setUDataAt(uDims[0] - 1, j, k, 0); // desna granica
        }
    }

    // Y pravac (V komponenta) - donja i gornja granica
    for (int i = 0; i < vDims[0]; i++)
    {
        for (int k = 0; k < vDims[2]; k++)
        {
            velocity->setVDataAt(i, 0, k, 0);          // donja granica
            velocity->setVDataAt(i, vDims[1] - 1, k, 0); // gornja granica
        }
    }

    // Z pravac (W komponenta) - prednja i zadnja granica
    for (int i = 0; i < wDims[0]; i++)
    {
        for (int j = 0; j < wDims[1]; j++)
        {
            velocity->setWDataAt(i, j, 0, 0);          // prednja granica
            velocity->setWDataAt(i, j, wDims[2] - 1, 0); // zadnja granica
        }
    }
}

void GridFluidSolver::applyColliderBoundaryConditions(std::shared_ptr<MACVectorGrid> velocity)
{
    Size3 uDims = velocity->getUDims();
    Size3 vDims = velocity->getVDims();
    Size3 wDims = velocity->getWDims();
    Size3 res = solverGrids.res();

    // ========================================================================
    // U KOMPONENTA (x-pravac)
    // ========================================================================
    for (int i = 0; i < uDims[0]; i++) {
        for (int j = 0; j < uDims[1]; j++) {
            for (int k = 0; k < uDims[2]; k++) {
                Vec3 pos = velocity->dataCenterPosU(i, j, k);
                double sdfValue = colliderSDF->sample(pos);

                if (sdfValue <= 0.0) {
                    // Indeks ćelije koja sadrži ovo lice (leva ili desna)
                    int cell_i = std::max(0, std::min(res[0]-1, i));
                    int cell_j = j;
                    int cell_k = k;
                    
                    Vec3 normal = colliderSDF->gradientAt(cell_i, cell_j, cell_k);
                    double len = normal.norm();
                    if (len > 1e-8) normal /= len;
                    
                    // Interpolacija V i W na poziciju U lica
                    double u_vel = velocity->u(i, j, k);
                    double v_vel = velocity->sampleVAt(pos);  // treba implementirati
                    double w_vel = velocity->sampleWAt(pos);  // treba implementirati
                    
                    // Projekcija brzine na normalu
                    double vn = u_vel * normal[0] + v_vel * normal[1] + w_vel * normal[2];
                    
                    // Oduzmi normalnu komponentu (free-slip)
                    double u_new = u_vel - vn * normal[0];
                    
                    velocity->setUDataAt(i, j, k, u_new); // free-slip
					//velocity->setUDataAt(i, j, k, 0);       // no-slip
                }
            }
        }
    }

    // ========================================================================
    // V KOMPONENTA (y-pravac) - slično
    // ========================================================================
    for (int i = 0; i < vDims[0]; i++) {
        for (int j = 0; j < vDims[1]; j++) {
            for (int k = 0; k < vDims[2]; k++) {
                Vec3 pos = velocity->dataCenterPosV(i, j, k);
                double sdfValue = colliderSDF->sample(pos);

                if (sdfValue <= 0.0) {
                    int cell_i = i;
                    int cell_j = std::max(0, std::min(res[1]-1, j));
                    int cell_k = k;
                    
                    Vec3 normal = colliderSDF->gradientAt(cell_i, cell_j, cell_k);
                    double len = normal.norm();
                    if (len > 1e-8) normal /= len;
                    
                    double u_vel = velocity->sampleUAt(pos);
                    double v_vel = velocity->v(i, j, k);
                    double w_vel = velocity->sampleWAt(pos);
                    
                    double vn = u_vel * normal[0] + v_vel * normal[1] + w_vel * normal[2];
                    double v_new = v_vel - vn * normal[1];
                    
					velocity->setVDataAt(i, j, k, v_new); // free-slip
					//velocity->setVDataAt(i, j, k, 0);       // no-slip
                }
            }
        }
    }

    // ========================================================================
    // W KOMPONENTA (z-pravac) - slično
    // ========================================================================
    for (int i = 0; i < wDims[0]; i++) {
        for (int j = 0; j < wDims[1]; j++) {
            for (int k = 0; k < wDims[2]; k++) {
                Vec3 pos = velocity->dataCenterPosW(i, j, k);
                double sdfValue = colliderSDF->sample(pos);

                if (sdfValue <= 0.0) {
                    int cell_i = i;
                    int cell_j = j;
                    int cell_k = std::max(0, std::min(res[2]-1, k));
                    
                    Vec3 normal = colliderSDF->gradientAt(cell_i, cell_j, cell_k);
                    double len = normal.norm();
                    if (len > 1e-8) normal /= len;
                    
                    double u_vel = velocity->sampleUAt(pos);
                    double v_vel = velocity->sampleVAt(pos);
                    double w_vel = velocity->w(i, j, k);
                    
                    double vn = u_vel * normal[0] + v_vel * normal[1] + w_vel * normal[2];
                    double w_new = w_vel - vn * normal[2];
                    
					velocity->setWDataAt(i, j, k, w_new);  // free-slip
					// velocity->setWDataAt(i, j, k, 0);       // no-slip
                }
            }
        }
    }

    // ========================================================================
    // 2. KOREKCIJA LEVEL SET-a (fluid ne sme biti unutar kolajdera)
    // ========================================================================
    
    for (int i = 0; i < res[0]; ++i) {
        for (int j = 0; j < res[1]; ++j) {
            for (int k = 0; k < res[2]; ++k) {
                Vec3 pos = solverGrids.getScalarGrid("FREE_SURFACE")->cellCenterPosition(Size3(i, j, k));
                double phi_collider = colliderSDF->sample(pos);
                
                if (phi_collider <= 0.0) solverGrids.getScalarGrid("FREE_SURFACE")->setDataAt(Size3(i, j, k), 1);
            }
        }
    }
}

Size3 GridFluidSolver::res() const
{ return solverRes; }
Vec3 GridFluidSolver::origin() const
{ return solverOrigin; }
Vec3 GridFluidSolver::spacing() const
{ return solverSpacing; }

GridManager &GridFluidSolver::grids()
{
    return solverGrids;
}

void GridFluidSolver::reinitializeFreeSurface(const ScalarGrid &inputMesh, double maxDistance, ScalarGrid &outputMesh)
{
    CenterScalarGrid tempGrid(inputMesh.res(), inputMesh.origin(), inputMesh.spacing(), 0);

    inputMesh.foreachIndex([&](int i, int j, int k)
                           {
                               outputMesh.setDataAt(Size3(i, j, k), inputMesh.at(i, j, k));
                               tempGrid.setDataAt(Size3(i, j, k), inputMesh.at(i, j, k));
                           });

    double pseudoT = Utils::pseudoTimeStep(inputMesh, inputMesh.spacing(), 0.5);
    int numIters = static_cast<int>(std::ceil(maxDistance / pseudoT));
    Vec3 spacing = inputMesh.spacing();
    double minSpacing = std::min(spacing[0], std::min(spacing[1], spacing[2]));
    for (int iter = 0; iter < numIters; iter++)
    {
        inputMesh.foreachIndex([&](int i, int j, int k)
                               {
                                   double sdf = outputMesh.at(i, j, k);
                                   double sign = sdf / std::sqrt(sdf * sdf + minSpacing * minSpacing);

                                   std::array<double, 2> minDx, minDy, minDz;

                                   //get derivatives
                                   getDerivatives(outputMesh, Size3(i, j, k), minDx, minDy, minDz);
                                   std::array<double, 2> maxDx = minDx, maxDy = minDy, maxDz = minDz;

                                   minDx[0] = std::max(0.0, minDx[0]);
                                   minDx[1] = std::min(0.0, minDx[1]);
                                   minDy[0] = std::max(0.0, minDy[0]);
                                   minDy[1] = std::min(0.0, minDy[1]);
                                   minDz[0] = std::max(0.0, minDz[0]);
                                   minDz[1] = std::min(0.0, minDz[1]);

                                   minDx[0] *= minDx[0];
                                   minDx[1] *= minDx[1];
                                   minDy[0] *= minDy[0];
                                   minDy[1] *= minDy[1];
                                   minDz[0] *= minDz[0];
                                   minDz[1] *= minDz[1];

                                   maxDx[0] = std::min(0.0, maxDx[0]);
                                   maxDx[1] = std::max(0.0, maxDx[1]);
                                   maxDy[0] = std::min(0.0, maxDy[0]);
                                   maxDy[1] = std::max(0.0, maxDy[1]);
                                   maxDz[0] = std::min(0.0, maxDz[0]);
                                   maxDz[1] = std::max(0.0, maxDz[1]);

                                   maxDx[0] *= maxDx[0];
                                   maxDx[1] *= maxDx[1];
                                   maxDy[0] *= maxDy[0];
                                   maxDy[1] *= maxDy[1];
                                   maxDz[0] *= maxDz[0];
                                   maxDz[1] *= maxDz[1];

                                   double minDistance =
                                       std::sqrt(
                                           minDx[0] + minDx[1] + minDy[0] + minDy[1] + minDz[0] + minDz[1]) - 1;
                                   double maxDistance =
                                       std::sqrt(
                                           maxDx[0] + maxDx[1] + maxDy[0] + maxDy[1] + maxDz[0] + maxDz[1]) - 1;

                                   double value = outputMesh.at(i,j,k)
                                       -pseudoT * std::max(sign, 0.0) * minDistance
                                       - pseudoT * std::min(sign, 0.0) * maxDistance;

                                   tempGrid.setDataAt(Size3(i, j, k), value);

                               });

        //swap grid values
        outputMesh.foreachIndex([&](int i, int j, int k)
                                {
                                    double tempValue = tempGrid.at(i, j, k);
                                    double value = outputMesh.at(i, j, k);

                                    outputMesh.setDataAt(Size3(i, j, k), tempValue);
                                    tempGrid.setDataAt(Size3(i, j, k), value);
                                });
    }
    std::cout << std::endl;
}

void GridFluidSolver::getDerivatives(ScalarGrid &grid,
                                     const Size3 indices,
                                     std::array<double, 2> &dx,
                                     std::array<double, 2> &dy,
                                     std::array<double, 2> &dz)
{
    auto upwindMethod = [](Vec3 values, double d) -> std::array<double,2>
    {
        double invd = 1 / d;
        std::array<double, 2> df{};
        df[0] = invd * (values[1] - values[0]);
        df[1] = invd * (values[2] - values[1]);
        return df;
    };

    int i = indices[0];
    int j = indices[1];
    int k = indices[2];

    Vec3 D0;
    Size3 size = grid.dataSize();
    Vec3 spacing = grid.spacing();

    const int im1 = (i < 1) ? 0 : i - 1;
    const int ip1 = std::min(i + 1, size[0] - 1);
    const int jm1 = (j < 1) ? 0 : j - 1;
    const int jp1 = std::min(j + 1, size[1] - 1);
    const int km1 = (k < 1) ? 0 : k - 1;
    const int kp1 = std::min(k + 1, size[2] - 1);

    D0[0] = grid.at(im1, j, k);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(ip1, j, k);
    dx = upwindMethod(D0, spacing[0]);

    D0[0] = grid.at(i, jm1, k);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(i, jp1, k);
    dy = upwindMethod(D0, spacing[1]);

    D0[0] = grid.at(i, j, km1);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(i, j, kp1);
    dz = upwindMethod(D0, spacing[2]);
}
void GridFluidSolver::reinitializeFreeSurface()
{
    std::cout << "reinitializeFreeSurface ..." << std::endl;
  
    std::shared_ptr<CenterScalarGrid>
        fluid = std::dynamic_pointer_cast<CenterScalarGrid>(grids().getScalarGrid("FREE_SURFACE"));
    assert(fluid != nullptr);

    CenterScalarGrid oldSdf(*fluid);

    double h = std::min(oldSdf.spacing()[0],std::min(oldSdf.spacing()[1],oldSdf.spacing()[2]));
    double minDistance  = 10.0;
    reinitializeFreeSurface(oldSdf, minDistance * h, *fluid);
}
void GridFluidSolver::extrapolateIntoCollider(ScalarGrid &grid)
{
    std::cout << "extrapolateIntoCollider ..." << std::endl;

    Size3 dims = grid.dataSize();
    unsigned int size = dims[0] * dims[1] * dims[2];

    //use to quickly determine which parts of the sdf grid is inside, with "true" meaning that the point is inside the sdf (a.k.a distance < 0)
    std::vector<bool> outsideMatrix(size);


    grid.foreachIndex([&](int i, int j, int k)
                      {
                          int flatIndex = k + dims[2] * (j + i * dims[1]);
                          //for now, evething is outside the boundaries, since there are no boundaries to worry about. update here later when you need
                          //actual boundaries
                          outsideMatrix[flatIndex] = true;
                      });

    std::vector<double> flatData;
    grid.getData(&flatData);
    std::vector<double> output(flatData.size());

    unsigned int depth = 1; //w/ CFL of 0.5 ceiling-ed up
    Utils::extrapolateToRegion(flatData, dims, outsideMatrix, depth, output);

    grid.foreachIndex([&](int i, int j, int k)
                      {
                          int flatIndex = k + dims[2] * (j + i * dims[1]);
                          grid.setDataAt(Size3(i, j, k), output[flatIndex]);
                      });
}

void GridFluidSolver::extrapolateVelocityToAir()
{
    std::cout << "extrapolateVelocityToAir ..." << std::endl;

    double cfl = 5.0;
    auto sdf = grids().getScalarGrid("FREE_SURFACE");
    auto vel = std::dynamic_pointer_cast<MACVectorGrid>(grids().getVelocityGrid());
    assert(vel != nullptr);

    std::vector<bool> insideMatrixU(static_cast<unsigned int>(vel->sizeU()));
    std::vector<bool> insideMatrixV(static_cast<unsigned int>(vel->sizeV()));
    std::vector<bool> insideMatrixW(static_cast<unsigned int>(vel->sizeW()));


    vel->forEachU([&](int i, int j, int k)
                  {
                      const Vec3 &position = vel->dataCenterPosU(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixU[vel->flatIndexU(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setUDataAt(i, j, k, 0);
                  });

    vel->forEachV([&](int i, int j, int k)
                  {
                      const Vec3 &position = vel->dataCenterPosV(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixV[vel->flatIndexV(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setVDataAt(i, j, k, 0);
                  });
    vel->forEachW([&](int i, int j, int k)
                  {
                      const Vec3 &position = vel->dataCenterPosW(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixW[vel->flatIndexW(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setWDataAt(i, j, k, 0);
                  });


    std::vector<double> uData, vData, wData;
    vel->getDataU(uData);
    vel->getDataV(vData);
    vel->getDataW(wData);

    std::vector<double> uOutput(uData.size()), vOutput(vData.size()), wOutput(wData.size());

    Utils::extrapolateToRegion(uData, vel->getUDims(), insideMatrixU, cfl, uOutput);
    Utils::extrapolateToRegion(vData, vel->getVDims(), insideMatrixV, cfl, vOutput);
    Utils::extrapolateToRegion(wData, vel->getWDims(), insideMatrixW, cfl, wOutput);

    vel->forEachU([&](int i, int j, int k)
                  {
                      double value = uOutput[vel->flatIndexU(i, j, k)];
                      vel->setUDataAt(i, j, k, value);
                  });

    vel->forEachV([&](int i, int j, int k)
                  {
                      double value = vOutput[vel->flatIndexV(i, j, k)];
                      vel->setVDataAt(i, j, k, value);
                  });
    vel->forEachW([&](int i, int j, int k)
                  {
                      double value = wOutput[vel->flatIndexW(i, j, k)];
                      vel->setWDataAt(i, j, k, value);
                  });

    applyBoundaryCondition();
}

void GridFluidSolver::resetSourceVelocity(double cylinderRadius,
    double cylinderHeight,
    double cylinderOffset,
    double sourceVelocity)
{
    // Dohvati potrebne grid-ove
    //std::cout << "resetSourceVelocity: getting velocity grid..." << std::endl;
    //std::cout.flush();
    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    
    //std::cout << "resetSourceVelocity: getting FREE_SURFACE grid..." << std::endl;
    std::cout.flush();
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("FREE_SURFACE"));
    
    //std::cout << "resetSourceVelocity: getting SOURCE_MASK grid..." << std::endl;
    std::cout.flush();
    auto sourceMask = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("SOURCE_MASK"));

    if (!velocity || !sdf)
    {
        printf("resetSourceVelocity: ERROR - velocity or SDF grid not found!\n");
        return;
    }

    if (!sourceMask)
    {
        printf("resetSourceVelocity: WARNING - SOURCE_MASK grid not found! Source region will not be protected.\n");
    }

    Size3 res = solverGrids.res();
    Vec3 spacing = solverGrids.spacing();
    double dx = spacing[0], dy = spacing[1], dz = spacing[2];

    double Lx = res[0] * dx;
    double Ly = res[1] * dy;
    double Lz = res[2] * dz;

    //std::cout << "resetSourceVelocity: Lx=" << Lx << " Ly=" << Ly << " Lz=" << Lz << std::endl;
    //std::cout.flush();

    // Centri cilindara
    Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);
    Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);

    //std::cout << "resetSourceVelocity: cyl1_center=(" << cyl1_center[0] << "," << cyl1_center[1] << "," << cyl1_center[2] << ")" << std::endl;
    //std::cout << "resetSourceVelocity: cyl2_center=(" << cyl2_center[0] << "," << cyl2_center[1] << "," << cyl2_center[2] << ")" << std::endl;
    //std::cout.flush();

    // ========================================================================
    // 1. POSTAVLJANJE BRZINE NA U LICIMA (X pravac - prvi cilindar)
    // ========================================================================
    Size3 uDims = velocity->getUDims();
    //std::cout << "resetSourceVelocity: uDims = (" << uDims[0] << "," << uDims[1] << "," << uDims[2] << ")" << std::endl;
    //std::cout.flush();

    int setU = 0;
    int sampleU = 0;

    for (int i = 0; i < uDims[0]; ++i)
    {
        for (int j = 0; j < uDims[1]; ++j)
        {
            for (int k = 0; k < uDims[2]; ++k)
            {
                Vec3 pos = velocity->dataCenterPosU(i, j, k);
                double x = pos[0], y = pos[1], z = pos[2];

                double dy1 = y - cyl1_center[1];
                double dz1 = z - cyl1_center[2];
                double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
                double axial_dist1 = x - cylinderOffset;
                bool inside1 = (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight);

                if (inside1)
                {
                    velocity->setUDataAt(i, j, k, sourceVelocity);
                    ++setU;

                    if (Scenario == 6 && sampleU < 4)
                    {
                        //printf("resetSourceVelocity U sample: i=%d j=%d k=%d pos=(%.3f,%.3f,%.3f) val=%.3f\n",
                        //       i, j, k, x, y, z, sourceVelocity);
                        ++sampleU;
                    }
                }
            }
        }
    }
    //std::cout << "resetSourceVelocity: U loop done, setU=" << setU << std::endl;
    //std::cout.flush();

    // ========================================================================
    // 2. POSTAVLJANJE BRZINE NA W LICIMA (Z pravac - drugi cilindar)
    // ========================================================================
    Size3 wDims = velocity->getWDims();
    //std::cout << "resetSourceVelocity: wDims = (" << wDims[0] << "," << wDims[1] << "," << wDims[2] << ")" << std::endl;
    //std::cout.flush();

    int setW = 0;
    int sampleW = 0;

    for (int i = 0; i < wDims[0]; ++i)
    {
        for (int j = 0; j < wDims[1]; ++j)
        {
            for (int k = 0; k < wDims[2]; ++k)
            {
                Vec3 pos = velocity->dataCenterPosW(i, j, k);
                double x = pos[0], y = pos[1], z = pos[2];

                double dx2 = x - cyl2_center[0];
                double dy2 = y - cyl2_center[1];
                double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
                double axial_dist2 = (Lz - cylinderOffset) - z;
                bool inside2 = (radial_dist2 < cylinderRadius && axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight);

                if (inside2)
                {
                    velocity->setWDataAt(i, j, k, -sourceVelocity);
                    ++setW;

                    if (Scenario == 6 && sampleW < 4)
                    {
                        //printf("resetSourceVelocity W sample: i=%d j=%d k=%d pos=(%.3f,%.3f,%.3f) val=%.3f\n", i, j, k, x, y, z, -sourceVelocity);
                        ++sampleW;
                    }
                }
            }
        }
    }
    //std::cout << "resetSourceVelocity: W loop done, setW=" << setW << std::endl;
    //std::cout.flush();

// ========================================================================
// 3. POSTAVLJANJE SDF-a (FREE_SURFACE) I SOURCE_MASK-a
// ========================================================================
//std::cout << "resetSourceVelocity: starting SDF loop, res=(" << res[0] << "," << res[1] << "," << res[2] << ")" << std::endl;
std::cout.flush();

//std::cout << "resetSourceVelocity: checking sourceMask pointer = " << sourceMask.get() << std::endl;
//std::cout.flush();

int sdfMarked = 0;
int sourceMaskMarked = 0;

//std::cout << "resetSourceVelocity: entering first for loop (i)" << std::endl;
//std::cout.flush();

for (int i = 0; i < res[0]; ++i)
{
    if (i % 10 == 0) {
        //std::cout << "  i=" << i << std::endl;
        //std::cout.flush();
    }
    
    for (int j = 0; j < res[1]; ++j)
    {
        for (int k = 0; k < res[2]; ++k)
        {
            // Prva iteracija - poseban ispis
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "  First cell (0,0,0): calculating position..." << std::endl;
                //std::cout.flush();
            }
            
            double cx = (i + 0.5) * dx;
            double cy = (j + 0.5) * dy;
            double cz = (k + 0.5) * dz;
            
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "    cx=" << cx << " cy=" << cy << " cz=" << cz << std::endl;
                //std::cout.flush();
            }
            
            // Provera za prvi cilindar
            double dy1 = cy - cyl1_center[1];
            double dz1 = cz - cyl1_center[2];
            double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
            double axial_dist1 = cx - cylinderOffset;
            bool inside1 = (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight);
            
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "    dy1=" << dy1 << " dz1=" << dz1 << " radial_dist1=" << radial_dist1 << " axial_dist1=" << axial_dist1 << " inside1=" << inside1 << std::endl;
                //std::cout.flush();
            }
            
            // Provera za drugi cilindar
            double dx2 = cx - cyl2_center[0];
            double dy2 = cy - cyl2_center[1];
            double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
            double axial_dist2 = (Lz - cylinderOffset) - cz;
            bool inside2 = (radial_dist2 < cylinderRadius && axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight);
            
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "    dx2=" << dx2 << " dy2=" << dy2 << " radial_dist2=" << radial_dist2 << " axial_dist2=" << axial_dist2 << " inside2=" << inside2 << std::endl;
                //std::cout.flush();
            }
            
            bool isInsideSource = inside1 || inside2;
            
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "    isInsideSource=" << isInsideSource << std::endl;
                //std::cout.flush();
                //std::cout << "    First cell - before setting values" << std::endl;
                //std::cout.flush();
            }

            if (isInsideSource)
            {
                if (i == 0 && j == 0 && k == 0)
                {
                    //std::cout << "    Calling sdf->setDataAt..." << std::endl;
                    //std::cout.flush();
                }
                
                sdf->setDataAt(Size3(i, j, k), -1.0);
                sdfMarked++;

                if (sourceMask)
                {
                    if (i == 0 && j == 0 && k == 0)
                    {
                        //std::cout << "    Calling sourceMask->setDataAt (inside)..." << std::endl;
                        //std::cout.flush();
                    }
                    sourceMask->setDataAt(Size3(i, j, k), -1.0);
                    sourceMaskMarked++;
                }
            }
            else
            {
                if (sourceMask)
                {
                    if (i == 0 && j == 0 && k == 0)
                    {
                        //std::cout << "    Setting sourceMask to 1.0 for (0,0,0)" << std::endl;
                        //std::cout.flush();
                    }
                    sourceMask->setDataAt(Size3(i, j, k), 1.0);

                    if (i == 0 && j == 0 && k == 0)
                    {
                        //std::cout << "    sourceMask->setDataAt done" << std::endl;
                        //std::cout.flush();
                    }
                }
            }
            
            if (i == 0 && j == 0 && k == 0)
            {
                //std::cout << "    First cell DONE" << std::endl;
                //std::cout.flush();
            }
        }
    }
}
#if 0
std::cout << "resetSourceVelocity: SDF loop done" << std::endl;
std::cout.flush();

    // ========================================================================
    // 4. DEBUG ISPIS ZA PRVIH NEKOLIKO ĆELIJA UNUTAR CILINDRA
    // ========================================================================
    if (Scenario == 6)
    {
        int checkCount = 0;
        printf("\n=== DEBUG: Prvih 5 celija unutar cilindra 1 ===\n");
        for (int i = 0; i < res[0] && checkCount < 5; ++i)
        {
            for (int j = 0; j < res[1] && checkCount < 5; ++j)
            {
                for (int k = 0; k < res[2] && checkCount < 5; ++k)
                {
                    double cx = (i + 0.5) * dx;
                    double cy = (j + 0.5) * dy;
                    double cz = (k + 0.5) * dz;
                    
                    double dy1 = cy - cyl1_center[1];
                    double dz1 = cz - cyl1_center[2];
                    double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
                    double axial_dist1 = cx - cylinderOffset;
                    bool inside1 = (radial_dist1 < cylinderRadius && axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight);
                    
                    if (inside1)
                    {
                        double sdfValue = sdf->sample(Vec3(cx, cy, cz));
                        double maskValue = sourceMask ? sourceMask->sample(Vec3(cx, cy, cz)) : 999;
                        printf("  Cell [%d,%d,%d] pos=(%.3f,%.3f,%.3f): radial=%.3f, axial=%.3f, SDF=%.3f, sourceMask=%.3f\n",
                               i, j, k, cx, cy, cz, radial_dist1, axial_dist1, sdfValue, maskValue);
                        checkCount++;
                    }
                }
            }
        }
        
        if (checkCount == 0)
        {
            printf("ERROR: Nijedna celija nije detektovana unutar cilindra 1!\n");
            printf("  Cyl1 center: (%.3f, %.3f, %.3f)\n", cyl1_center[0], cyl1_center[1], cyl1_center[2]);
            printf("  Radius: %.3f\n", cylinderRadius);
        }
    }

    // ========================================================================
    // 5. ZAVRŠNI DEBUG ISPISI
    // ========================================================================
    if (Scenario == 6)
    {
        printf("resetSourceVelocity: setU=%d setW=%d\n", setU, setW);
        printf("resetSourceVelocity: SDF marked %d cells as fluid\n", sdfMarked);
        if (sourceMask)
        {
            printf("resetSourceVelocity: SOURCE_MASK marked %d cells\n", sourceMaskMarked);
        }
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << "resetSourceVelocity: END - SUCCESS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout.flush();
#endif
}
