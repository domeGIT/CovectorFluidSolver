#include "SemiLagrangianAdvectSolver.h"
#include <cmath>
#include <iostream>
#include <algorithm>

// ============================================================================
// POMONE FUNKCIJE
// ============================================================================

double safe_val(double x)
{
    if (std::isnan(x) || std::isinf(x)) return 0.0;
    if (x > 1e6) return 1e6;
    if (x < -1e6) return -1e6;
    return x;
}

Vec3 safe_val(const Vec3& v)
{
    return Vec3(safe_val(v[0]), safe_val(v[1]), safe_val(v[2]));
}

Mat3::Mat3()
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m[i][j] = (i == j ? 1.0 : 0.0);
}

Vec3 Mat3::operator*(const Vec3& v) const
{
    return Vec3(
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
    );
}

Mat3 mat3_transpose(const Mat3& A)
{
    Mat3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T.m[i][j] = A.m[j][i];
    return T;
}

// ============================================================================
// KLASA IMPLEMENTACIJA
// ============================================================================

SemiLagrangianAdvectSolver::SemiLagrangianAdvectSolver() {}
SemiLagrangianAdvectSolver::~SemiLagrangianAdvectSolver() {}

// ----------------------------------------------------------------------------
// RK4 Backtrace (Algorithm 3, Line 2)
// ----------------------------------------------------------------------------
Vec3 SemiLagrangianAdvectSolver::backtrace_rk4(
    std::shared_ptr<MACVectorGrid> flow,
    double dt,
    Vec3 start,
    std::shared_ptr<ScalarGrid> colliderSDF)
{
    Vec3 pt = start;
    double remaining = dt;
    Vec3 spc = flow->spacing();
    double h = std::min(spc[0], std::min(spc[1], spc[2]));

    while (remaining > 1e-8)
    {
        double step = std::min(remaining, h * 0.5);

        // RK4 integration (unazad)
        Vec3 k1 = flow->sample(pt);
        Vec3 k2 = flow->sample(pt - 0.5 * step * k1);
        Vec3 k3 = flow->sample(pt - 0.5 * step * k2);
        Vec3 k4 = flow->sample(pt - step * k3);

        Vec3 avgVel = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        Vec3 newPt = pt - step * avgVel;

        // Boundary handling (Section 5.4.5)
        double phi0 = colliderSDF->sample(pt);
        double phi1 = colliderSDF->sample(newPt);

        if (phi0 < 0.0) {
            double w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1) + 1.e-12);
            newPt = w * pt + (1.0 - w) * newPt;
            pt = newPt;
            break;
        }

        pt = newPt;
        remaining -= step;
    }

    return pt;
}

// ----------------------------------------------------------------------------
// Sample MAC Velocity (trilinear interpolation)
// ----------------------------------------------------------------------------
Vec3 SemiLagrangianAdvectSolver::sample_mac_velocity(
    std::shared_ptr<MACVectorGrid> grid,
    double x, double y, double z)
{
    auto lerp = [](double a, double b, double t) { return (1.0 - t) * a + t * b; };

    Vec3 spc = grid->spacing();
    double dx = spc[0], dy = spc[1], dz = spc[2];
    Size3 dim = grid->res();
    int Nx = dim[0], Ny = dim[1], Nz = dim[2];

    // Clamp na domen (FACE-based granice)
    x = std::clamp(x, 0.0, Nx * dx);
    y = std::clamp(y, 0.0, Ny * dy);
    z = std::clamp(z, 0.0, Nz * dz);

    // U-komponenta (x-faces)
    int iu = static_cast<int>(std::floor(x / dx));
    int ju = static_cast<int>(std::floor((y - 0.5 * dy) / dy));
    int ku = static_cast<int>(std::floor((z - 0.5 * dz) / dz));
    iu = std::clamp(iu, 0, Nx - 1);
    ju = std::clamp(ju, 0, Ny - 1);
    ku = std::clamp(ku, 0, Nz - 1);
    double fu = std::clamp((x / dx) - iu, 0.0, 1.0);
    double u_val = lerp(grid->u(iu, ju, ku), grid->u(iu + 1, ju, ku), fu);

    // V-komponenta (y-faces)
    int iv = static_cast<int>(std::floor((x - 0.5 * dx) / dx));
    int jv = static_cast<int>(std::floor(y / dy));
    int kv = static_cast<int>(std::floor((z - 0.5 * dz) / dz));
    iv = std::clamp(iv, 0, Nx - 1);
    jv = std::clamp(jv, 0, Ny - 1);
    kv = std::clamp(kv, 0, Nz - 1);
    double fv = std::clamp((y / dy) - jv, 0.0, 1.0);
    double v_val = lerp(grid->v(iv, jv, kv), grid->v(iv, jv + 1, kv), fv);

    // W-komponenta (z-faces)
    int iw = static_cast<int>(std::floor((x - 0.5 * dx) / dx));
    int jw = static_cast<int>(std::floor((y - 0.5 * dy) / dy));
    int kw = static_cast<int>(std::floor(z / dz));
    iw = std::clamp(iw, 0, Nx - 1);
    jw = std::clamp(jw, 0, Ny - 1);
    kw = std::clamp(kw, 0, Nz - 1);
    double fw = std::clamp((z / dz) - kw, 0.0, 1.0);
    double w_val = lerp(grid->w(iw, jw, kw), grid->w(iw, jw, kw + 1), fw);

    return Vec3(u_val, v_val, w_val);
}

// ----------------------------------------------------------------------------
// MAC-Aware Jacobian Transpose (Eq. 40)
// ----------------------------------------------------------------------------
Mat3 SemiLagrangianAdvectSolver::compute_jacobian_transpose_MAC_aware(
    std::shared_ptr<MACVectorGrid> flow,
    double dt,
    double x, double y, double z,
    std::shared_ptr<ScalarGrid> colliderSDF,
    bool& success,
    int faceAxis)
{
    success = false;

    Vec3 spc = flow->spacing();
    double dx = spc[0], dy = spc[1], dz = spc[2];

    // PROVERA 1: Da li smo u fluidu?
    double phi = colliderSDF->sample(Vec3(x, y, z));
    if (phi > dx * 0.5) {
        success = true;
        return Mat3();  // Identity za air
    }

    // PROVERA 2: Blizu slobodne povrine? (Section 5.4.5)
    if (phi > -2.0 * dx && phi < 2.0 * dx) {
        success = true;
        return Mat3();  // Fallback na SL
    }

    Vec3 dPsi_dx, dPsi_dy, dPsi_dz;

    // ========================================================================
    // MAC-AWARE FINITE DIFFERENCES prema Eq. (40)
    // Normalni izvod: cell-centered differences
    // Tangencijalni izvodi: face-centered perturbacije
    // ========================================================================

    if (faceAxis == 0) {  // U-FACE: normalni pravac je X
        // /x: cell centers na x  0.5*dx
        Vec3 psi_plus = this->backtrace_rk4(flow, dt, Vec3(x + 0.5 * dx, y, z), colliderSDF);
        Vec3 psi_minus = this->backtrace_rk4(flow, dt, Vec3(x - 0.5 * dx, y, z), colliderSDF);
        dPsi_dx = (psi_plus - psi_minus) / dx;

        // /y, /z: face perturbacije
        double h = dy * 0.5;
        Vec3 psi_yp = this->backtrace_rk4(flow, dt, Vec3(x, y + h, z), colliderSDF);
        Vec3 psi_ym = this->backtrace_rk4(flow, dt, Vec3(x, y - h, z), colliderSDF);
        dPsi_dy = (psi_yp - psi_ym) / (2.0 * h);

        h = dz * 0.5;
        Vec3 psi_zp = this->backtrace_rk4(flow, dt, Vec3(x, y, z + h), colliderSDF);
        Vec3 psi_zm = this->backtrace_rk4(flow, dt, Vec3(x, y, z - h), colliderSDF);
        dPsi_dz = (psi_zp - psi_zm) / (2.0 * h);
    }
    else if (faceAxis == 1) {  // V-FACE: normalni pravac je Y
        // /y: cell centers
        Vec3 psi_plus = this->backtrace_rk4(flow, dt, Vec3(x, y + 0.5 * dy, z), colliderSDF);
        Vec3 psi_minus = this->backtrace_rk4(flow, dt, Vec3(x, y - 0.5 * dy, z), colliderSDF);
        dPsi_dy = (psi_plus - psi_minus) / dy;

        // /x, /z: face perturbacije
        double h = dx * 0.5;
        Vec3 psi_xp = this->backtrace_rk4(flow, dt, Vec3(x + h, y, z), colliderSDF);
        Vec3 psi_xm = this->backtrace_rk4(flow, dt, Vec3(x - h, y, z), colliderSDF);
        dPsi_dx = (psi_xp - psi_xm) / (2.0 * h);

        h = dz * 0.5;
        Vec3 psi_zp = this->backtrace_rk4(flow, dt, Vec3(x, y, z + h), colliderSDF);
        Vec3 psi_zm = this->backtrace_rk4(flow, dt, Vec3(x, y, z - h), colliderSDF);
        dPsi_dz = (psi_zp - psi_zm) / (2.0 * h);
    }
    else {  // W-FACE: normalni pravac je Z
        // /z: cell centers
        Vec3 psi_plus = this->backtrace_rk4(flow, dt, Vec3(x, y, z + 0.5 * dz), colliderSDF);
        Vec3 psi_minus = this->backtrace_rk4(flow, dt, Vec3(x, y, z - 0.5 * dz), colliderSDF);
        dPsi_dz = (psi_plus - psi_minus) / dz;

        // /x, /y: face perturbacije
        double h = dx * 0.5;
        Vec3 psi_xp = this->backtrace_rk4(flow, dt, Vec3(x + h, y, z), colliderSDF);
        Vec3 psi_xm = this->backtrace_rk4(flow, dt, Vec3(x - h, y, z), colliderSDF);
        dPsi_dx = (psi_xp - psi_xm) / (2.0 * h);

        h = dy * 0.5;
        Vec3 psi_yp = this->backtrace_rk4(flow, dt, Vec3(x, y + h, z), colliderSDF);
        Vec3 psi_ym = this->backtrace_rk4(flow, dt, Vec3(x, y - h, z), colliderSDF);
        dPsi_dy = (psi_yp - psi_ym) / (2.0 * h);
    }

    // Konstruii Jacobian D
    Mat3 D_psi;
    D_psi.m[0][0] = dPsi_dx[0]; D_psi.m[0][1] = dPsi_dy[0]; D_psi.m[0][2] = dPsi_dz[0];
    D_psi.m[1][0] = dPsi_dx[1]; D_psi.m[1][1] = dPsi_dy[1]; D_psi.m[1][2] = dPsi_dz[1];
    D_psi.m[2][0] = dPsi_dx[2]; D_psi.m[2][1] = dPsi_dy[2]; D_psi.m[2][2] = dPsi_dz[2];

    // PROVERA 3: Determinanta (Section 5.4.3)
    double det = D_psi.m[0][0] * (D_psi.m[1][1] * D_psi.m[2][2] - D_psi.m[1][2] * D_psi.m[2][1]) -
        D_psi.m[0][1] * (D_psi.m[1][0] * D_psi.m[2][2] - D_psi.m[1][2] * D_psi.m[2][0]) +
        D_psi.m[0][2] * (D_psi.m[1][0] * D_psi.m[2][1] - D_psi.m[1][1] * D_psi.m[2][0]);

    if (std::abs(det) < 1e-8) {
        success = true;
        std::cout << "Fallback!" << std::endl;
        return Mat3();  // Fallback na identity
    }

    success = true;
    return mat3_transpose(D_psi);
}

// ----------------------------------------------------------------------------
// BFECC Limiter (Section 5.4.2)
// ----------------------------------------------------------------------------
Vec3 SemiLagrangianAdvectSolver::bfecc_limiter_MAC(
    const Vec3& value,
    double x, double y, double z,
    std::shared_ptr<MACVectorGrid> velocity_ref,
    std::shared_ptr<ScalarGrid> colliderSDF,
    double dx, double dy, double dz,
    int component)
{
    // Nai cell-centered indekse
    int ci = static_cast<int>(std::floor(x / dx));
    int cj = static_cast<int>(std::floor(y / dy));
    int ck = static_cast<int>(std::floor(z / dz));

    Size3 dim = velocity_ref->res();
    int Nx = dim[0], Ny = dim[1], Nz = dim[2];

    // Nai min/max u 3x3x3 susedstvu
    double minVal = value[component];
    double maxVal = value[component];

    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            for (int dk = -1; dk <= 1; ++dk) {
                int i = ci + di;
                int j = cj + dj;
                int k = ck + dk;

                if (i >= 0 && i < Nx && j >= 0 && j < Ny && k >= 0 && k < Nz) {
                    double nx = (i + 0.5) * dx;
                    double ny = (j + 0.5) * dy;
                    double nz = (k + 0.5) * dz;

                    double phi = colliderSDF->sample(Vec3(nx, ny, nz));
                    if (phi < 0.0) {  // Samo fluid elije
                        Vec3 neighbor = this->sample_mac_velocity(velocity_ref, nx, ny, nz);
                        minVal = std::min(minVal, neighbor[component]);
                        maxVal = std::max(maxVal, neighbor[component]);
                    }
                }
            }
        }
    }

    // Clamp na min/max (bez irenja opsega!)
    Vec3 result = value;
    result[component] = std::clamp(value[component], minVal, maxVal);

    return result;
}

// ----------------------------------------------------------------------------
// 1-Step Covector Advection (Algorithm 3)
// ----------------------------------------------------------------------------
void SemiLagrangianAdvectSolver::covector_advection_1step(
    std::shared_ptr<MACVectorGrid> velocity_in,
    double dt,
    std::shared_ptr<MACVectorGrid> velocity_out,
    std::shared_ptr<MACVectorGrid> flow,
    std::shared_ptr<ScalarGrid> colliderSDF)
{
    Size3 dim = velocity_in->res();
    int Nx = dim[0], Ny = dim[1], Nz = dim[2];
    Vec3 spc = velocity_in->spacing();
    double dx = spc[0], dy = spc[1], dz = spc[2];

    // U-komponenta (x-faces)
#pragma omp parallel for
    for (int i = 0; i < Nx + 1; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                double x = i * dx;
                double y = (j + 0.5) * dy;
                double z = (k + 0.5) * dz;

                // 1. Backtrace do (x)
                Vec3 psi = this->backtrace_rk4(flow, dt, Vec3(x, y, z), colliderSDF);

                // 2. Interpoliraj ceo vektor na (x)
                Vec3 u_psi = this->sample_mac_velocity(velocity_in, psi[0], psi[1], psi[2]);

                // 3. Izraunaj (D) sa MAC-aware Jacobian
                bool jac_success = false;
                Mat3 D_psi_T = this->compute_jacobian_transpose_MAC_aware(
                    flow, dt, x, y, z, colliderSDF, jac_success, 0);

                if (!jac_success) {
                    velocity_out->setUDataAt(i, j, k, velocity_in->u(i, j, k));
                    continue;
                }

                // 4. Primeni transformaciju: u_new = (D)  u()
                Vec3 u_corrected = D_psi_T * u_psi;
                velocity_out->setUDataAt(i, j, k, safe_val(u_corrected[0]));
            }
        }
    }

    // V-komponenta (y-faces)
#pragma omp parallel for
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny + 1; ++j) {
            for (int k = 0; k < Nz; ++k) {
                double x = (i + 0.5) * dx;
                double y = j * dy;
                double z = (k + 0.5) * dz;

                Vec3 psi = this->backtrace_rk4(flow, dt, Vec3(x, y, z), colliderSDF);
                Vec3 u_psi = this->sample_mac_velocity(velocity_in, psi[0], psi[1], psi[2]);

                bool jac_success = false;
                Mat3 D_psi_T = this->compute_jacobian_transpose_MAC_aware(
                    flow, dt, x, y, z, colliderSDF, jac_success, 1);

                if (!jac_success) {
                    velocity_out->setVDataAt(i, j, k, velocity_in->v(i, j, k));
                    continue;
                }

                Vec3 v_corrected = D_psi_T * u_psi;
                velocity_out->setVDataAt(i, j, k, safe_val(v_corrected[1]));
            }
        }
    }

    // W-komponenta (z-faces)
#pragma omp parallel for
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz + 1; ++k) {
                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;
                double z = k * dz;

                Vec3 psi = this->backtrace_rk4(flow, dt, Vec3(x, y, z), colliderSDF);
                Vec3 u_psi = this->sample_mac_velocity(velocity_in, psi[0], psi[1], psi[2]);

                bool jac_success = false;
                Mat3 D_psi_T = this->compute_jacobian_transpose_MAC_aware(
                    flow, dt, x, y, z, colliderSDF, jac_success, 2);

                if (!jac_success) {
                    velocity_out->setWDataAt(i, j, k, velocity_in->w(i, j, k));
                    continue;
                }

                Vec3 w_corrected = D_psi_T * u_psi;
                velocity_out->setWDataAt(i, j, k, safe_val(w_corrected[2]));
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Covector Correction sa BFECC (Algorithm 4)
// ----------------------------------------------------------------------------
void SemiLagrangianAdvectSolver::covector_correction_BFECC(
    std::shared_ptr<MACVectorGrid> velocity_prev,
    double dt,
    std::shared_ptr<MACVectorGrid> velocity_SL,
    std::shared_ptr<MACVectorGrid> velocity,
    std::shared_ptr<ScalarGrid> colliderSDF)
{
    std::cout << "===== COVECTOR CORRECTION (BFECC) =====" << std::endl;

    // ========================================================================
    // Korak 1: Forward covector advection (Algorithm 4, Line 1)
    // ========================================================================
    std::cout << "  Step 1: Forward covector advection..." << std::endl;
    this->covector_advection_1step(velocity_prev, dt, velocity, velocity_prev, colliderSDF);

    // ========================================================================
    // Korak 2: Backward advection za error estimation (Algorithm 4, Line 2)
    // ========================================================================
    std::cout << "  Step 2: Backward advection..." << std::endl;
    auto velocity_temp = std::make_shared<MACVectorGrid>(*velocity);
    auto velocity_back = std::make_shared<MACVectorGrid>(*velocity_prev);
    this->covector_advection_1step(velocity, -dt, velocity_back, velocity_prev, colliderSDF);

    // ========================================================================
    // Korak 3: Izraunaj roundtrip error (Algorithm 4, Line 3)
    // ========================================================================
    std::cout << "  Step 3: Computing roundtrip error..." << std::endl;
    auto error_grid = std::make_shared<MACVectorGrid>(*velocity_prev);
    error_grid->forEachU([&](int i, int j, int k) {
        double err = velocity_back->u(i, j, k) - velocity_prev->u(i, j, k);
        error_grid->setUDataAt(i, j, k, err);
        });
    error_grid->forEachV([&](int i, int j, int k) {
        double err = velocity_back->v(i, j, k) - velocity_prev->v(i, j, k);
        error_grid->setVDataAt(i, j, k, err);
        });
    error_grid->forEachW([&](int i, int j, int k) {
        double err = velocity_back->w(i, j, k) - velocity_prev->w(i, j, k);
        error_grid->setWDataAt(i, j, k, err);
        });

    // ========================================================================
    // Korak 4: Error correction (Algorithm 4, Line 4)
    // ========================================================================
    std::cout << "  Step 4: Error correction..." << std::endl;
    auto error_corrected = std::make_shared<MACVectorGrid>(*error_grid);
    this->covector_advection_1step(error_grid, dt, error_corrected, velocity_prev, colliderSDF);

    // Oduzmi korigovanu greku od velocity (e/2)
    velocity->forEachU([&](int i, int j, int k) {
        double u_new = velocity->u(i, j, k) - 0.5 * error_corrected->u(i, j, k);
        velocity->setUDataAt(i, j, k, safe_val(u_new));
        });
    velocity->forEachV([&](int i, int j, int k) {
        double v_new = velocity->v(i, j, k) - 0.5 * error_corrected->v(i, j, k);
        velocity->setVDataAt(i, j, k, safe_val(v_new));
        });
    velocity->forEachW([&](int i, int j, int k) {
        double w_new = velocity->w(i, j, k) - 0.5 * error_corrected->w(i, j, k);
        velocity->setWDataAt(i, j, k, safe_val(w_new));
        });

    std::cout << "===== COVECTOR CORRECTION COMPLETE =====" << std::endl;
}

// ----------------------------------------------------------------------------
// Advect (MACVectorGrid) - Standardna SL advekcija
// ----------------------------------------------------------------------------
void SemiLagrangianAdvectSolver::advect(MACVectorGrid& input,
    std::shared_ptr<MACVectorGrid> flow,
    double dt,
    MACVectorGrid& output,
    std::shared_ptr<ScalarGrid> colliderSDF)
{
    std::cout << "advect velocity SL..." << std::endl;

    double h = std::min(input.spacing()[0], std::min(input.spacing()[1], input.spacing()[2]));

    input.forEachU([&](int i, int j, int k) {
        Vec3 inPos = input.dataCenterPosU(i, j, k);
        Vec3 outPos = output.dataCenterPosU(i, j, k);
        if (colliderSDF->sample(inPos) > 0.0) {
            Vec3 outPoint = this->backtrace_rk4(flow, dt, outPos, colliderSDF);
            Vec3 sample = input.sample(outPoint);
            output.setUDataAt(i, j, k, sample[0]);
        }
        });

    input.forEachV([&](int i, int j, int k) {
        Vec3 inPos = input.dataCenterPosV(i, j, k);
        Vec3 outPos = output.dataCenterPosV(i, j, k);
        if (colliderSDF->sample(inPos) > 0.0) {
            Vec3 outPoint = this->backtrace_rk4(flow, dt, outPos, colliderSDF);
            Vec3 sample = input.sample(outPoint);
            output.setVDataAt(i, j, k, sample[1]);
        }
        });

    input.forEachW([&](int i, int j, int k) {
        Vec3 inPos = input.dataCenterPosW(i, j, k);
        Vec3 outPos = output.dataCenterPosW(i, j, k);
        if (colliderSDF->sample(inPos) > 0.0) {
            Vec3 outPoint = this->backtrace_rk4(flow, dt, outPos, colliderSDF);
            output.setWDataAt(i, j, k, input.sample(outPoint)[2]);
        }
        });
}

// ----------------------------------------------------------------------------
// Advect (ScalarGrid) - Standardna SL advekcija
// ----------------------------------------------------------------------------
void SemiLagrangianAdvectSolver::advect(ScalarGrid& input,
    std::shared_ptr<MACVectorGrid> flow,
    double dt,
    ScalarGrid& output,
    std::shared_ptr<ScalarGrid> colliderSDF)
{
    std::cout << "advect surface SL ..." << std::endl;

    double h = std::min(output.spacing()[0], std::min(output.spacing()[1], output.spacing()[2]));
    Size3 size = input.dataSize();

#pragma omp parallel for
    for (int k = 0; k < (int)size[2]; k++)
        for (int j = 0; j < (int)size[1]; j++)
            for (int i = 0; i < (int)size[0]; i++) {
                Vec3 inPos = input.dataCenterPosition(i, j, k);
                double collideSample = colliderSDF->sample(inPos);
                collideSample = std::isnan(collideSample) ? 1.0 : collideSample;

                if (collideSample > 0.0) {
                    Vec3 outPoint = this->backtrace_rk4(flow, dt, inPos, colliderSDF);
                    double outPointSample = input.cubicSample(outPoint);
                    output.setDataAt(Size3(i, j, k), outPointSample);
                }
            }
}

double SemiLagrangianAdvectSolver::compute_kinetic_energy(
    std::shared_ptr<MACVectorGrid> velocity)
{
    Size3 dim = velocity->res();
    int Nx = dim(0), Ny = dim(1), Nz = dim(2);
    Vec3 spc = velocity->spacing();
    double dx = spc(0), dy = spc(1), dz = spc(2);
    double cell_volume = dx * dy * dz;

    double kinetic_energy = 0.0;

    // Za svaku eliju (cell-centered)
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // Pozicija elijskog centra
                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;
                double z = (k + 0.5) * dz;

                // Interpoliraj brzinu u centru elije
                Vec3 v = this->sample_mac_velocity(velocity, x, y, z);

                // Kinetika energija = 0.5 * |v|^2 * zapremina
                double cell_energy = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * cell_volume;
                kinetic_energy += cell_energy;
            }
        }
    }

    return kinetic_energy;
}

Vec3 SemiLagrangianAdvectSolver::compute_curl_at_cell(
    std::shared_ptr<MACVectorGrid> velocity,
    int i, int j, int k)  // cell-centered indeksi
{
    Size3 dim = velocity->res();
    int Nx = dim(0), Ny = dim(1), Nz = dim(2);
    Vec3 spc = velocity->spacing();
    double dx = spc(0), dy = spc(1), dz = spc(2);

    // Ogranii indekse za centralnu diferenciju
    int im1 = std::max(i - 1, 0);
    int ip1 = std::min(i + 1, Nx - 1);
    int jm1 = std::max(j - 1, 0);
    int jp1 = std::min(j + 1, Ny - 1);
    int km1 = std::max(k - 1, 0);
    int kp1 = std::min(k + 1, Nz - 1);

    // Izraunaj sve tri komponente vorticity-ja:
    // _x = w/y - v/z
    // _y = u/z - w/x
    // _z = v/x - u/y

    // Za _x trebaju nam v i w na y i z face-ovima
    double dw_dy, dv_dz;

    // w/y - koristimo w na z-face-ovima
    double w_yp = velocity->w(i, jp1, k) + velocity->w(i, jp1, k + 1); // interpolacija?
    double w_ym = velocity->w(i, jm1, k) + velocity->w(i, jm1, k + 1);
    w_yp *= 0.5; w_ym *= 0.5;  // proseno na elijski centar
    dw_dy = (w_yp - w_ym) / (2.0 * dy);

    // v/z - koristimo v na y-face-ovima
    double v_zp = velocity->v(i, j, kp1) + velocity->v(i + 1, j, kp1); // interpolacija
    double v_zm = velocity->v(i, j, km1) + velocity->v(i + 1, j, km1);
    v_zp *= 0.5; v_zm *= 0.5;
    dv_dz = (v_zp - v_zm) / (2.0 * dz);

    double wx = dw_dy - dv_dz;

    // _y = u/z - w/x
    double du_dz, dw_dx;

    // u/z
    double u_zp = velocity->u(ip1, j, kp1) + velocity->u(ip1, j + 1, kp1); // interpolacija
    double u_zm = velocity->u(ip1, j, km1) + velocity->u(ip1, j + 1, km1);
    u_zp *= 0.5; u_zm *= 0.5;
    du_dz = (u_zp - u_zm) / (2.0 * dz);

    // w/x
    double w_xp = velocity->w(ip1, j, k) + velocity->w(ip1, j, k + 1);
    double w_xm = velocity->w(im1, j, k) + velocity->w(im1, j, k + 1);
    w_xp *= 0.5; w_xm *= 0.5;
    dw_dx = (w_xp - w_xm) / (2.0 * dx);

    double wy = du_dz - dw_dx;

    // _z = v/x - u/y
    double dv_dx, du_dy;

    // v/x
    double v_xp = velocity->v(ip1, j, k) + velocity->v(ip1, j + 1, k); // interpolacija
    double v_xm = velocity->v(im1, j, k) + velocity->v(im1, j + 1, k);
    v_xp *= 0.5; v_xm *= 0.5;
    dv_dx = (v_xp - v_xm) / (2.0 * dx);

    // u/y
    double u_yp = velocity->u(i, jp1, k) + velocity->u(i + 1, jp1, k);
    double u_ym = velocity->u(i, jm1, k) + velocity->u(i + 1, jm1, k);
    u_yp *= 0.5; u_ym *= 0.5;
    du_dy = (u_yp - u_ym) / (2.0 * dy);

    double wz = dv_dx - du_dy;

    return Vec3(wx, wy, wz);
}

double SemiLagrangianAdvectSolver::compute_enstrophy(
    std::shared_ptr<MACVectorGrid> velocity)
{
    Size3 dim = velocity->res();
    int Nx = dim(0), Ny = dim(1), Nz = dim(2);
    Vec3 spc = velocity->spacing();
    double dx = spc(0), dy = spc(1), dz = spc(2);
    double cell_volume = dx * dy * dz;

    double enstrophy = 0.0;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                Vec3 curl = compute_curl_at_cell(velocity, i, j, k);

                // Enstrophy = ||^2 * zapremina
                double cell_enstrophy = (curl[0] * curl[0] +
                    curl[1] * curl[1] +
                    curl[2] * curl[2]) * cell_volume;
                enstrophy += cell_enstrophy;
            }
        }
    }

    return enstrophy;
}

double SemiLagrangianAdvectSolver::compute_enstrophy_simple(
    std::shared_ptr<MACVectorGrid> velocity)
{
    Size3 dim = velocity->res();
    int Nx = dim(0), Ny = dim(1), Nz = dim(2);
    Vec3 spc = velocity->spacing();
    double dx = spc(0), dy = spc(1), dz = spc(2);
    double cell_volume = dx * dy * dz;

    double enstrophy = 0.0;

    for (int k = 1; k < Nz - 1; ++k) {
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                // Jednostavnija centralna diferencija
                double dw_dy = (velocity->w(i, j + 1, k) - velocity->w(i, j - 1, k)) / (2 * dy);
                double dv_dz = (velocity->v(i, j, k + 1) - velocity->v(i, j, k - 1)) / (2 * dz);
                double wx = dw_dy - dv_dz;

                double du_dz = (velocity->u(i + 1, j, k + 1) - velocity->u(i + 1, j, k - 1)) / (2 * dz);
                double dw_dx = (velocity->w(i + 1, j, k) - velocity->w(i - 1, j, k)) / (2 * dx);
                double wy = du_dz - dw_dx;

                double dv_dx = (velocity->v(i + 1, j, k) - velocity->v(i - 1, j, k)) / (2 * dx);
                double du_dy = (velocity->u(i, j + 1, k) - velocity->u(i, j - 1, k)) / (2 * dy);
                double wz = dv_dx - du_dy;

                double vort2 = wx * wx + wy * wy + wz * wz;
                enstrophy += vort2 * cell_volume;
            }
        }
    }

    return enstrophy;
}

bool SemiLagrangianAdvectSolver::is_degenerate_jacobian(const Mat3& J)
{
    // Proveri da li je Jacobian previe blizu nule ili ima ekstremne vrednosti
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (std::isnan(J.m[i][j]) || std::isinf(J.m[i][j]) || std::abs(J.m[i][j]) > 10.0) {
                return true;
            }
        }
    }
    return false;
}

