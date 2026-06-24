#ifndef SEMILAGRANGIANADVECTSOLVER_H
#define SEMILAGRANGIANADVECTSOLVER_H

#include "MACVectorGrid.h"
#include "ScalarGrid.h"
#include <memory>
#include <vector>
#include <limits>
#include <algorithm>

struct Mat3
{
    double m[3][3];
    Mat3();
    Vec3 operator*(const Vec3& v) const;
};

Mat3 mat3_transpose(const Mat3& A);

class SemiLagrangianAdvectSolver
{
public:
    SemiLagrangianAdvectSolver();
    virtual ~SemiLagrangianAdvectSolver();

    void advect(MACVectorGrid& input,
        std::shared_ptr<MACVectorGrid> flow,
        double dt,
        MACVectorGrid& output,
        std::shared_ptr<ScalarGrid> colliderSDF);

    void advect(ScalarGrid& input,
        std::shared_ptr<MACVectorGrid> flow,
        double dt,
        ScalarGrid& output,
        std::shared_ptr<ScalarGrid> colliderSDF);

    // Glavna funkcija za kovektorsku korekciju sa BFECC
    void covector_correction_BFECC(
        std::shared_ptr<MACVectorGrid> velocity_prev,
        double dt,
        std::shared_ptr<MACVectorGrid> velocity_SL,
        std::shared_ptr<MACVectorGrid> velocity,
        std::shared_ptr<ScalarGrid> colliderSDF);

    // Pomona: 1-step covector advection (Algorithm 3)
    void covector_advection_1step(
        std::shared_ptr<MACVectorGrid> velocity_in,
        double dt,
        std::shared_ptr<MACVectorGrid> velocity_out,
        std::shared_ptr<MACVectorGrid> flow,
        std::shared_ptr<ScalarGrid> colliderSDF);

    // Kinetic energy and enstrophy
    double SemiLagrangianAdvectSolver::compute_kinetic_energy(
        std::shared_ptr<MACVectorGrid> velocity);

    double SemiLagrangianAdvectSolver::compute_enstrophy(
        std::shared_ptr<MACVectorGrid> velocity);

    double SemiLagrangianAdvectSolver::compute_enstrophy_simple(
        std::shared_ptr<MACVectorGrid> velocity);



private:
    Vec3 backtrace_rk4(
        std::shared_ptr<MACVectorGrid> flow,
        double dt,
        Vec3 start,
        std::shared_ptr<ScalarGrid> colliderSDF);

    Vec3 sample_mac_velocity(
        std::shared_ptr<MACVectorGrid> grid,
        double x, double y, double z);

    Mat3 compute_jacobian_transpose_MAC_aware(
        std::shared_ptr<MACVectorGrid> flow,
        double dt,
        double x, double y, double z,
        std::shared_ptr<ScalarGrid> colliderSDF,
        bool& success,
        int faceAxis);  // 0=u-face, 1=v-face, 2=w-face

    Vec3 bfecc_limiter_MAC(
        const Vec3& value,
        double x, double y, double z,
        std::shared_ptr<MACVectorGrid> velocity_ref,
        std::shared_ptr<ScalarGrid> colliderSDF,
        double dx, double dy, double dz,
        int component);

    bool  is_degenerate_jacobian(const Mat3& J);

    Vec3 SemiLagrangianAdvectSolver::compute_curl_at_cell(
        std::shared_ptr<MACVectorGrid> velocity,
        int i, int j, int k);
};

#endif // SEMILAGRANGIANADVECTSOLVER_H