#include "PressureSolver.h"
#include<eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Array<int, 3, 1> Size3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;

PressureSolver::PressureSolver() {}
PressureSolver::~PressureSolver() {}

void PressureSolver::solve(const MACVectorGrid &input,
                           const ScalarGrid &colliderSDF,
                           const ScalarGrid &fluidSDF,
                           MACVectorGrid &output,
                           double dt,
                           const ScalarGrid* sourceMask)
{
    //std::cout << "PressureSolver::solve - START" << std::endl;

    setupMarkers(input, colliderSDF, fluidSDF, sourceMask);
    //std::cout << "PressureSolver::solve - after setupMarkers" << std::endl;

    setupSystem(input, dt);
    //std::cout << "PressureSolver::solve - after setupSystem" << std::endl;


    //Symmetric matrix, so the solver only has to resolve the upper part

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper > solver;
    x = solver.compute(A).solve(b);
    //std::cout << "PressureSolver::solve - after CG solve" << std::endl;

    applyPressureGradient(input, output, dt);
    //std::cout << "PressureSolver::solve - END" << std::endl;
}

void PressureSolver::setupMarkers(const MACVectorGrid &input,
                                  const ScalarGrid &collider,
                                  const ScalarGrid &fluid,
                                  const ScalarGrid* sourceMask)
{
    Size3 size = input.res();

    markers.clear();
    markers.resize(size[0] * size[1] * size[2]);
    sourceRegion.clear();
    sourceRegion.resize(size[0] * size[1] * size[2], 0);

    #pragma omp parallel for
    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {
                int flatIndex = k + j * size[2] + i * size[2] * size[1];

                Vec3 pos = input.cellCenterPosition(Size3(i, j, k));
                double fluidSample = fluid.sample(pos);
                
                if (collider.sample(pos) < 0.0)
                {
                    markers[flatIndex] = COLLIDER;
                }
                else if (sourceMask != nullptr && sourceMask->sample(pos) < 0.0)
                {
                    // mark source region as FLUID with special behavior
                    markers[flatIndex] = FLUID;
                    sourceRegion[flatIndex] = 1;
                }
                else if (fluidSample < 0.0)
                {
                    markers[flatIndex] = FLUID;
                }
                else
                {
                    markers[flatIndex] = AIR;
                }
            }
        }
    }

    // ========== DEBUG ISPIS - VAN PARALELNE PETLJE ==========
    int sourceRegionCount = 0;
    for (int idx = 0; idx < sourceRegion.size(); idx++)
    {
        if (sourceRegion[idx]) sourceRegionCount++;
    }
    //std::cout << "PressureSolver::setupMarkers: sourceRegion count = " << sourceRegionCount << std::endl;
}

void PressureSolver::setupSystem(const MACVectorGrid &input, double dt)
{
    typedef Eigen::Triplet<double> DataPoint;

    Size3 res = input.res();
    int size = res[0] * res[1] * res[2];
    A = Eigen::SparseMatrix<double>(size, size);
    std::vector<DataPoint> Adata;
    b = Eigen::VectorXd::Zero(size);
    x.resize(size);

    Vec3 c(1 / (input.spacing()[0] * input.spacing()[0]),
           1 / (input.spacing()[1] * input.spacing()[1]),
           1 / (input.spacing()[2] * input.spacing()[2]));

    int iIter = res[1] * res[2];
    int jIter = res[2];
    int kIter = 1;

    for (int diagIndex = 0; diagIndex < size; diagIndex++)
    {
        int unflatI = diagIndex / (res[1] * res[2]);
        int temp = diagIndex - unflatI * res[1] * res[2];
        int unflatJ = temp / res[2];
        int unflatK = temp % res[2];

        // SLUČAJ 1: Source region (inlet cilindar)
        if (markers[diagIndex] == FLUID && sourceRegion[diagIndex])
        {
            // Dirichlet: p = 0
            Adata.emplace_back(diagIndex, diagIndex, 1.0);
            b[diagIndex] = 0.0;
            continue;
        }
        
        // SLUČAJ 2: Normalan FLUID (nije source region)
        else if (markers[diagIndex] == FLUID)
        {
            double value = input.divergenceAt(Size3(unflatI, unflatJ, unflatK));
            b[diagIndex] = -value / dt;
            
            double diag_val = 0.0;
            
            // +X smer (desno)
            if (unflatI + 1 < res[0])
            {
                int neighborIdx = diagIndex + iIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[0];
                    // Dodaj off-diagonal samo ako sused NIJE source region
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[0]);
                    }
                }
            }
            
            // -X smer (levo)
            if (unflatI > 0)
            {
                int neighborIdx = diagIndex - iIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[0];
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[0]);
                    }
                }
            }
            
            // +Y smer (gore)
            if (unflatJ + 1 < res[1])
            {
                int neighborIdx = diagIndex + jIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[1];
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[1]);
                    }
                }
            }
            
            // -Y smer (dole)
            if (unflatJ > 0)
            {
                int neighborIdx = diagIndex - jIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[1];
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[1]);
                    }
                }
            }
            
            // +Z smer (napred)
            if (unflatK + 1 < res[2])
            {
                int neighborIdx = diagIndex + kIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[2];
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[2]);
                    }
                }
            }
            
            // -Z smer (nazad)
            if (unflatK > 0)
            {
                int neighborIdx = diagIndex - kIter;
                if (markers[neighborIdx] != COLLIDER)
                {
                    diag_val += c[2];
                    if (markers[neighborIdx] == FLUID && !sourceRegion[neighborIdx])
                    {
                        Adata.emplace_back(diagIndex, neighborIdx, -c[2]);
                    }
                }
            }
            
            Adata.emplace_back(diagIndex, diagIndex, diag_val);
        }
        else
        {
            // COLLIDER ili AIR
            Adata.emplace_back(diagIndex, diagIndex, 1.0);
            b[diagIndex] = 0.0;
        }
    }

    A.setFromTriplets(Adata.begin(), Adata.end());
}

void PressureSolver::applyPressureGradient(const MACVectorGrid& input,
    MACVectorGrid& output,
    double dt)
{
    Size3 res = input.res();
    Vec3 invSpacing(1.0 / input.spacing()[0],
        1.0 / input.spacing()[1],
        1.0 / input.spacing()[2]);

    // Update U faces (i from 0 to res[0]) -> U dims = res[0]+1, res[1], res[2]

    Size3 uDims = input.getUDims();
    #pragma omp parallel for
    for (int i = 0; i < uDims[0]; ++i)
    {
        for (int j = 0; j < uDims[1]; ++j)
        {
            for (int k = 0; k < uDims[2]; ++k)
            {
                // u-face at (i,j,k) is between cells (i-1,j,k) and (i,j,k) in cell-index space
                int leftI = i - 1;
                int rightI = i;
                bool leftValid = (leftI >= 0 && leftI < res[0]);
                bool rightValid = (rightI >= 0 && rightI < res[0]);

                int leftFlat = leftValid ? (k + j * res[2] + leftI * res[2] * res[1]) : -1;
                int rightFlat = rightValid ? (k + j * res[2] + rightI * res[2] * res[1]) : -1;

                // skip faces adjacent to collider cells
                if ((leftValid && markers[leftFlat] == COLLIDER) || (rightValid && markers[rightFlat] == COLLIDER))
                    continue;
                // If either adjacent cell is part of the source region, preserve input velocity (Dirichlet inflow)
                if ((leftValid && sourceRegion[leftFlat]) || (rightValid && sourceRegion[rightFlat]))
                {
                    output.setUDataAt(i, j, k, input.u(i, j, k));
                    continue;
                }
                double p_left = (leftValid) ? x[leftFlat] : x[rightFlat];
                double p_right = (rightValid) ? x[rightFlat] : x[leftFlat];
                double uData = input.u(i, j, k);
                double grad_p = (p_right - p_left) * invSpacing[0];
                output.setUDataAt(i, j, k, uData - dt * grad_p);
            }
        }
    }

    // Update V faces

    Size3 vDims = input.getVDims();
    #pragma omp parallel for
    for (int i = 0; i < vDims[0]; ++i)
    {
        for (int j = 0; j < vDims[1]; ++j)
        {
            for (int k = 0; k < vDims[2]; ++k)
            {
                int bottomJ = j - 1;
                int topJ = j;
                bool bottomValid = (bottomJ >= 0 && bottomJ < res[1]);
                bool topValid = (topJ >= 0 && topJ < res[1]);

                int bottomFlat = bottomValid ? (k + bottomJ * res[2] + i * res[2] * res[1]) : -1;
                int topFlat = topValid ? (k + topJ * res[2] + i * res[2] * res[1]) : -1;

                if ((bottomValid && markers[bottomFlat] == COLLIDER) || (topValid && markers[topFlat] == COLLIDER))
                    continue;
                if ((bottomValid && sourceRegion[bottomFlat]) || (topValid && sourceRegion[topFlat]))
                {
                    output.setVDataAt(i, j, k, input.v(i, j, k));
                    continue;
                }
                double p_bottom = (bottomValid) ? x[bottomFlat] : x[topFlat];
                double p_top = (topValid) ? x[topFlat] : x[bottomFlat];
                double vData = input.v(i, j, k);
                double grad_p = (p_top - p_bottom) * invSpacing[1];
                output.setVDataAt(i, j, k, vData - dt * grad_p);
            }
        }
    }

    // Update W faces

    Size3 wDims = input.getWDims();
    #pragma omp parallel for
    for (int i = 0; i < wDims[0]; ++i)
    {
        for (int j = 0; j < wDims[1]; ++j)
        {
            for (int k = 0; k < wDims[2]; ++k)
            {
                int backK = k - 1;
                int frontK = k;
                bool backValid = (backK >= 0 && backK < res[2]);
                bool frontValid = (frontK >= 0 && frontK < res[2]);

                int backFlat = backValid ? (backK + j * res[2] + i * res[2] * res[1]) : -1;
                int frontFlat = frontValid ? (frontK + j * res[2] + i * res[2] * res[1]) : -1;

                if ((backValid && markers[backFlat] == COLLIDER) || (frontValid && markers[frontFlat] == COLLIDER))
                    continue;
                if ((backValid && sourceRegion[backFlat]) || (frontValid && sourceRegion[frontFlat]))
                {
                    output.setWDataAt(i, j, k, input.w(i, j, k));
                    continue;
                }
                double p_back = (backValid) ? x[backFlat] : x[frontFlat];
                double p_front = (frontValid) ? x[frontFlat] : x[backFlat];
                double wData = input.w(i, j, k);
                double grad_p = (p_front - p_back) * invSpacing[2];
                output.setWDataAt(i, j, k, wData - dt * grad_p);
            }
        }
    }
}



