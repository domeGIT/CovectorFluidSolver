#include "MACVectorGrid.h"
#include "Utils.h"

extern int Scenario;
extern double sourceVelocity;
extern double dx;
extern double cylinderRadius;
extern double cylinderHeight;
extern double cylinderOffset;
extern double sourceVelocity;

Vec3 MACVectorGrid::sample(const Vec3 point)
{
// U MACVectorGrid::sample(Vec3 pos) ili u funkciji gde se uzorkuje brzina za advekciju
if (Scenario == 6)
{
    // Parametri cilindra (moraju biti dostupni ili globalni ili prosleđeni)
    double dx = spacing()[0];
    double dy = spacing()[1];
    double dz = spacing()[2];
    
    Size3 res = this->res();
    double Lx = res[0] * dx;
    double Ly = res[1] * dy;
    double Lz = res[2] * dz;
        
    // Centri cilindara (isti kao u resetSourceVelocity)
    Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);
    Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);
    
    // Tolerancija (10% ćelije)
    double tolerance = dx * 0.1;
    
    // ============================================================
    // CILINDAR 1 (X pravac)
    // ============================================================
    double dy1 = point[1] - cyl1_center[1];
    double dz1 = point[2] - cyl1_center[2];
    double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
    double axial_dist1 = point[0] - cylinderOffset;
    
    bool inside1 = (radial_dist1 < cylinderRadius + tolerance && 
                    axial_dist1 >= -tolerance && 
                    axial_dist1 <= cylinderHeight + tolerance);
    
    if (inside1)
    {
        return Vec3(sourceVelocity, 0.0, 0.0);
    }
    
    // ============================================================
    // CILINDAR 2 (Z pravac)
    // ============================================================
    double dx2 = point[0] - cyl2_center[0];
    double dy2 = point[1] - cyl2_center[1];
    double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
    double axial_dist2 = (Lz - cylinderOffset) - point[2];
    
    bool inside2 = (radial_dist2 < cylinderRadius + tolerance && 
                    axial_dist2 >= -tolerance && 
                    axial_dist2 <= cylinderHeight + tolerance);
    
    if (inside2)
    {
        return Vec3(0.0, 0.0, -sourceVelocity);
    }
}

   return Utils::trilinearInterpolateFromGrid(*this,point);
}

Vec3 MACVectorGrid::cubicSample(const Vec3 point)
{
    return Utils::cubicInterpolateFromGrid(*this,point);
}

double MACVectorGrid::divergenceAt(int i, int j, int k) const
{
    double left, right, front, back, up, down;

    left = u(i, j, k);
    right = u(i + 1, j, k);

    down = v(i, j, k);
    up = v(i, j + 1, k);

    back = w(i, j, k);
    front = w(i, j, k + 1);

    double x, y, z;

    x = (right - left) / spacing()[0];
    y = (up - down) / spacing()[1];
    z = (front - back) / spacing()[2];


    return x + y + z;
}

double MACVectorGrid::divergenceAt(Size3 indices) const
{
    return divergenceAt(indices[0], indices[1], indices[2]);
}

double MACVectorGrid::divergence(const Vec3 x) const
{
    double result = 0.0;
    std::array<Size3, 8> indices;
    std::array<double, 8> weights;

    Utils::getCoordinatesAndWeights(res() - Size3(1, 1, 1),
                                    spacing(),
                                    origin() + 0.5 * spacing(),
                                    x,
                                    indices,
                                    weights,
                                    nullptr);

    for (int i = 0; i < 8; i++)
    {
        result += weights[i] * divergenceAt(indices[i]);
    }

    return result;

}

Vec3 MACVectorGrid::curlAt(Size3 indices) const
{
    return curlAt(indices[0], indices[1], indices[2]);
}

Vec3 MACVectorGrid::curlAt(int i, int j, int k) const
{
    Vec3 left, right, front, back, up, down;
    const Size3 gridRes = res();
    const Vec3 gridSpacing = spacing();

    left = atCenter(std::max(i - 1, (int) 0), j, k);
    right = atCenter(std::min(i + 1, gridRes[0] - 1), j, k);

    up = atCenter(i, std::max(j - 1, (int) 0), k);
    down = atCenter(i, std::min(j + 1, gridRes[1] - 1), k);

    back = atCenter(i, j, std::max(k - 1, (int) 0));
    front = atCenter(i, j, std::min(k + 1, gridRes[2] - 1));

    double x, y, z;

    x = 0.5 * (up[2] - down[2]) / gridSpacing[1] - 0.5 * (front[1] - back[1]) / gridSpacing[2];
    y = 0.5 * (front[0] - back[0]) / gridSpacing[2] - 0.5 * (left[2] - right[2]) / gridSpacing[0];
    z = 0.5 * (left[1] - right[1]) / gridSpacing[0] - 0.5 * (up[0] - down[0]) / gridSpacing[1];

    return Vec3(x, y, z);
}

Vec3 MACVectorGrid::curl(const Vec3 point) const
{
    Vec3 result(0, 0, 0);
    std::array<Size3, 8> indices;
    std::array<double, 8> weights;

    Utils::getCoordinatesAndWeights(res() - Size3(1, 1, 1),
                                    spacing(),
                                    origin() + 0.5 * spacing(),
                                    point,
                                    indices,
                                    weights,
                                    nullptr);

    for (int i = 0; i < 8; i++)
    {
        result += weights[i] * curlAt(indices[i]);
    }

    return result;
}

Vec3 MACVectorGrid::atCenter(int i, int j, int k) const
{
    const Vec3 sample = Vec3(u(i, j, k) + u(i + 1, j, k), v(i, j, k) + v(i, j + 1, k),
                             w(i, j, k) + w(i, j, k + 1));
    return 0.5 * sample;
}

void MACVectorGrid::onResize(const Size3 &resolution, const Eigen::Vector3d gridSpacing, const Eigen::Vector3d origin,
                             Eigen::Vector3d initialValue)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    unsigned int sizeX = static_cast<unsigned int>(getUDims()[0] * getUDims()[1] * getUDims()[2]);
    unsigned int sizeY = static_cast<unsigned int>(getVDims()[0] * getVDims()[1] * getVDims()[2]);
    unsigned int sizeZ = static_cast<unsigned int>(getWDims()[0] * getWDims()[1] * getWDims()[2]);

    dataU.resize(sizeX);
    dataV.resize(sizeY);
    dataW.resize(sizeZ);

    std::fill(dataU.begin(),dataU.end(),initialValue[0]);
    std::fill(dataV.begin(),dataV.end(),initialValue[1]);
    std::fill(dataW.begin(),dataW.end(),initialValue[2]);

    double x = gridSpacing[0];
    double y = gridSpacing[1];
    double z = gridSpacing[2];

    originU = origin + 0.5 * Vec3(0, y, z);
    originV = origin + 0.5 * Vec3(x, 0, z);
    originW = origin + 0.5 * Vec3(x, y, 0);
}

double MACVectorGrid::divergence(int i, int j, int k) const
{
    return divergence(Vec3(i, j, k));
}

void MACVectorGrid::fillData(std::function<Vec3(Vec3)> predicate)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                Vec3 positionU = dataCenterPosU(i, j, k);
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(predicate(positionU)[0]);

                Vec3 positionV = dataCenterPosV(i, j, k);
                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(predicate(positionV)[1]);

                Vec3 positionW = dataCenterPosW(i, j, k);
                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(predicate(positionW)[2]);

            }
        }
    }
}

void MACVectorGrid::fillData(std::function<Vec3(double, double, double)> predicate)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                Vec3 positionU = dataCenterPosU(i, j, k);
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(predicate(positionU[0], positionU[1], positionU[2])[0]);

                Vec3 positionV = dataCenterPosV(i, j, k);
                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(predicate(positionV[0], positionV[1], positionV[2])[1]);

                Vec3 positionW = dataCenterPosW(i, j, k);
                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(predicate(positionW[0], positionW[1], positionW[2])[2]);

            }
        }
    }
}

void MACVectorGrid::fillData(double value)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(value);

                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(value);

                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(value);

            }
        }
    }
}

void MACVectorGrid::getDataU(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataU)
        result.push_back(value);
}

void MACVectorGrid::getDataV(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataV)
        result.push_back(value);
}

void MACVectorGrid::getDataW(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataW)
        result.push_back(value);
}

Vec3 MACVectorGrid::dataCenterPosW(const Vec3 indices) const
{
    return originW + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosW(int i, int j, int k) const
{
    return dataCenterPosW(Vec3(i, j, k));
}

Vec3 MACVectorGrid::dataCenterPosV(const Vec3 indices) const
{
    return originV + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosV(int i, int j, int k) const
{
    return dataCenterPosV(Vec3(i, j, k));
}

Vec3 MACVectorGrid::dataCenterPosU(const Vec3 indices) const
{
    return originU + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosU(int i, int j, int k) const
{
    return dataCenterPosU(Vec3(i, j, k));
}
double MACVectorGrid::sizeU()
{
    return dataU.size();
}

double MACVectorGrid::sizeV()
{
    return dataV.size();
}

double MACVectorGrid::sizeW()
{
    return dataW.size();
}

const double MACVectorGrid::u(int i, int j, int k) const
{
    return dataU[flatIndexU(i,j,k)];
}

const double MACVectorGrid::v(int i, int j, int k) const
{
    return dataV[flatIndexV(i,j,k)];
}

const double MACVectorGrid::w(int i, int j, int k) const
{
    return dataW[flatIndexW(i,j,k)];
}

double MACVectorGrid::u(int i, int j, int k)
{
    return dataU[flatIndexU(i,j,k)];
}

double MACVectorGrid::v(int i, int j, int k)
{
    return dataV[flatIndexV(i,j,k)];
}

double MACVectorGrid::w(int i, int j, int k)
{
    return dataW[flatIndexW(i,j,k)];
}
void MACVectorGrid::setUDataAt(int i, int j, int k, double value)
{
    int flatIdx = flatIndexU(i, j, k);
    dataU[flatIdx] = value;
}

void MACVectorGrid::setVDataAt(int i, int j, int k, double value)
{
    int flatIdx = flatIndexV(i, j, k);
    dataV[flatIdx] = value;
}

void MACVectorGrid::setWDataAt(int i, int j, int k, double value)
{
    int flatIdx = flatIndexW(i, j, k);;
    dataW[flatIdx] = value;
}

void MACVectorGrid::setUDataAt(Size3 indices, double value)
{
    setUDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::setVDataAt(Size3 indices, double value)
{
    setVDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::setWDataAt(Size3 indices, double value)
{
    setWDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::forEachU(std::function<void(int, int, int)> operation)
{
    #pragma omp parallel for
    for (int i = 0; i < getUDims()[0]; i++)
    {
        for (int j = 0; j < getUDims()[1]; j++)
        {
            for (int k = 0; k < getUDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}

void MACVectorGrid::forEachV(std::function<void(int, int, int)> operation)
{
    #pragma omp parallel for
    for (int i = 0; i < getVDims()[0]; i++)
    {
        for (int j = 0; j < getVDims()[1]; j++)
        {
            for (int k = 0; k < getVDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}

void MACVectorGrid::forEachW(std::function<void(int, int, int)> operation)
{
    #pragma omp parallel for
    for (int i = 0; i < getWDims()[0]; i++)
    {
        for (int j = 0; j < getWDims()[1]; j++)
        {
            for (int k = 0; k < getWDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}

double MACVectorGrid::maxNorm() const
{
    double maxNorm = 0.0;
    Size3 dim = this->res();  // ili dataSize(), zavisno od implementacije
    int Nx = dim[0], Ny = dim[1], Nz = dim[2];

    // u-komponenta (Nx+1, Ny, Nz)
    #pragma omp parallel for
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx + 1; ++i) {
                double u = this->u(i, j, k);
                maxNorm = std::max(maxNorm, std::abs(u));
            }

    // v-komponenta (Nx, Ny+1, Nz)
    #pragma omp parallel for
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny + 1; ++j)
            for (int i = 0; i < Nx; ++i) {
                double v = this->v(i, j, k);
                maxNorm = std::max(maxNorm, std::abs(v));
            }

    // w-komponenta (Nx, Ny, Nz+1)
    #pragma omp parallel for
    for (int k = 0; k < Nz + 1; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                double w = this->w(i, j, k);
                maxNorm = std::max(maxNorm, std::abs(w));
            }

    // Vrati max od svih komponenti (jo konzervativnije)
    return maxNorm;
}

// ============================================================================
// DODATNE FUNKCIJE ZA INTERPOLACIJU BRZINE NA PROIZVOLJNOJ POZICIJI
// ============================================================================

double MACVectorGrid::sampleUAt(const Vec3 point) const
{
    // Interpolacija U komponente (definisana na licima u X pravcu)
    // U komponenta je definisana na pozicijama (i*dx, (j+0.5)*dy, (k+0.5)*dz)
    
    Size3 res = this->res();
    Vec3 spacing = this->spacing();
    Vec3 origin = this->origin();
    
    double x = point[0];
    double y = point[1];
    double z = point[2];
    
    double dx = spacing[0];
    double dy = spacing[1];
    double dz = spacing[2];
    
    // Pronađi indekse okolnih U lica
    int i0 = std::max(0, std::min(res[0], (int)((x - origin[0]) / dx)));
    int i1 = std::min(i0 + 1, res[0]);
    
    int j0 = std::max(0, std::min(res[1] - 1, (int)((y - origin[1] - 0.5*dy) / dy)));
    int j1 = std::min(j0 + 1, res[1] - 1);
    
    int k0 = std::max(0, std::min(res[2] - 1, (int)((z - origin[2] - 0.5*dz) / dz)));
    int k1 = std::min(k0 + 1, res[2] - 1);
    
    // Uzorkovanje na 8 ćelija
    double u000 = u(i0, j0, k0);
    double u001 = u(i0, j0, k1);
    double u010 = u(i0, j1, k0);
    double u011 = u(i0, j1, k1);
    double u100 = u(i1, j0, k0);
    double u101 = u(i1, j0, k1);
    double u110 = u(i1, j1, k0);
    double u111 = u(i1, j1, k1);
    
    // Izračunaj težine za trilinearnu interpolaciju
    double x0 = origin[0] + i0 * dx;
    double x1 = origin[0] + i1 * dx;
    double y0 = origin[1] + (j0 + 0.5) * dy;
    double y1 = origin[1] + (j1 + 0.5) * dy;
    double z0 = origin[2] + (k0 + 0.5) * dz;
    double z1 = origin[2] + (k1 + 0.5) * dz;
    
    double tx = (x - x0) / (x1 - x0);
    double ty = (y - y0) / (y1 - y0);
    double tz = (z - z0) / (z1 - z0);
    
    // Trilinearna interpolacija
    double u00 = u000 * (1 - tx) + u100 * tx;
    double u01 = u001 * (1 - tx) + u101 * tx;
    double u10 = u010 * (1 - tx) + u110 * tx;
    double u11 = u011 * (1 - tx) + u111 * tx;
    
    double u0 = u00 * (1 - ty) + u10 * ty;
    double u1 = u01 * (1 - ty) + u11 * ty;
    
    return u0 * (1 - tz) + u1 * tz;
}

double MACVectorGrid::sampleVAt(const Vec3 point) const
{
    // Interpolacija V komponente (definisana na licima u Y pravcu)
    // V komponenta je definisana na pozicijama ((i+0.5)*dx, j*dy, (k+0.5)*dz)
    
    Size3 res = this->res();
    Vec3 spacing = this->spacing();
    Vec3 origin = this->origin();
    
    double x = point[0];
    double y = point[1];
    double z = point[2];
    
    double dx = spacing[0];
    double dy = spacing[1];
    double dz = spacing[2];
    
    // Pronađi indekse okolnih V lica
    int i0 = std::max(0, std::min(res[0] - 1, (int)((x - origin[0] - 0.5*dx) / dx)));
    int i1 = std::min(i0 + 1, res[0] - 1);
    
    int j0 = std::max(0, std::min(res[1], (int)((y - origin[1]) / dy)));
    int j1 = std::min(j0 + 1, res[1]);
    
    int k0 = std::max(0, std::min(res[2] - 1, (int)((z - origin[2] - 0.5*dz) / dz)));
    int k1 = std::min(k0 + 1, res[2] - 1);
    
    // Uzorkovanje na 8 ćelija
    double v000 = v(i0, j0, k0);
    double v001 = v(i0, j0, k1);
    double v010 = v(i0, j1, k0);
    double v011 = v(i0, j1, k1);
    double v100 = v(i1, j0, k0);
    double v101 = v(i1, j0, k1);
    double v110 = v(i1, j1, k0);
    double v111 = v(i1, j1, k1);
    
    // Izračunaj težine za trilinearnu interpolaciju
    double x0 = origin[0] + (i0 + 0.5) * dx;
    double x1 = origin[0] + (i1 + 0.5) * dx;
    double y0 = origin[1] + j0 * dy;
    double y1 = origin[1] + j1 * dy;
    double z0 = origin[2] + (k0 + 0.5) * dz;
    double z1 = origin[2] + (k1 + 0.5) * dz;
    
    double tx = (x - x0) / (x1 - x0);
    double ty = (y - y0) / (y1 - y0);
    double tz = (z - z0) / (z1 - z0);
    
    // Trilinearna interpolacija
    double v00 = v000 * (1 - tx) + v100 * tx;
    double v01 = v001 * (1 - tx) + v101 * tx;
    double v10 = v010 * (1 - tx) + v110 * tx;
    double v11 = v011 * (1 - tx) + v111 * tx;
    
    double v0 = v00 * (1 - ty) + v10 * ty;
    double v1 = v01 * (1 - ty) + v11 * ty;
    
    return v0 * (1 - tz) + v1 * tz;
}

double MACVectorGrid::sampleWAt(const Vec3 point) const
{
    // Interpolacija W komponente (definisana na licima u Z pravcu)
    // W komponenta je definisana na pozicijama ((i+0.5)*dx, (j+0.5)*dy, k*dz)
    
    Size3 res = this->res();
    Vec3 spacing = this->spacing();
    Vec3 origin = this->origin();
    
    double x = point[0];
    double y = point[1];
    double z = point[2];
    
    double dx = spacing[0];
    double dy = spacing[1];
    double dz = spacing[2];
    
    // Pronađi indekse okolnih W lica
    int i0 = std::max(0, std::min(res[0] - 1, (int)((x - origin[0] - 0.5*dx) / dx)));
    int i1 = std::min(i0 + 1, res[0] - 1);
    
    int j0 = std::max(0, std::min(res[1] - 1, (int)((y - origin[1] - 0.5*dy) / dy)));
    int j1 = std::min(j0 + 1, res[1] - 1);
    
    int k0 = std::max(0, std::min(res[2], (int)((z - origin[2]) / dz)));
    int k1 = std::min(k0 + 1, res[2]);
    
    // Uzorkovanje na 8 ćelija
    double w000 = w(i0, j0, k0);
    double w001 = w(i0, j0, k1);
    double w010 = w(i0, j1, k0);
    double w011 = w(i0, j1, k1);
    double w100 = w(i1, j0, k0);
    double w101 = w(i1, j0, k1);
    double w110 = w(i1, j1, k0);
    double w111 = w(i1, j1, k1);
    
    // Izračunaj težine za trilinearnu interpolaciju
    double x0 = origin[0] + (i0 + 0.5) * dx;
    double x1 = origin[0] + (i1 + 0.5) * dx;
    double y0 = origin[1] + (j0 + 0.5) * dy;
    double y1 = origin[1] + (j1 + 0.5) * dy;
    double z0 = origin[2] + k0 * dz;
    double z1 = origin[2] + k1 * dz;
    
    double tx = (x - x0) / (x1 - x0);
    double ty = (y - y0) / (y1 - y0);
    double tz = (z - z0) / (z1 - z0);
    
    // Trilinearna interpolacija
    double w00 = w000 * (1 - tx) + w100 * tx;
    double w01 = w001 * (1 - tx) + w101 * tx;
    double w10 = w010 * (1 - tx) + w110 * tx;
    double w11 = w011 * (1 - tx) + w111 * tx;
    
    double w0 = w00 * (1 - ty) + w10 * ty;
    double w1 = w01 * (1 - ty) + w11 * ty;
    
    return w0 * (1 - tz) + w1 * tz;
}

// ============================================================================
// INTERPOLACIJA CELOKUPNE BRZINE (U, V, W) NA PROIZVOLJNOJ POZICIJI
// ============================================================================

Vec3 MACVectorGrid::sampleVelocityAt(const Vec3 point) const
{
    return Vec3(sampleUAt(point), sampleVAt(point), sampleWAt(point));
}

#if 0
// ============================================================================
// 
// ============================================================================

Vec3 MACVectorGrid::vorticityVector(int i, int j, int k) const
{
    const Size3 gridRes = res();
    const Vec3 gridSpacing = spacing();
    double dx = gridSpacing[0];
    double dy = gridSpacing[1];
    double dz = gridSpacing[2];

    // Helper       
    auto safe_u = [&](int ii, int jj, int kk) -> double {
        ii = std::clamp(ii, 0, (int)getUDims()[0] - 1);
        return u(ii, jj, kk);
        };
    auto safe_v = [&](int ii, int jj, int kk) -> double {
        jj = std::clamp(jj, 0, (int)getVDims()[1] - 1);
        return v(ii, jj, kk);
        };
    auto safe_w = [&](int ii, int jj, int kk) -> double {
        kk = std::clamp(kk, 0, (int)getWDims()[2] - 1);
        return w(ii, jj, kk);
        };

    //       (i+0.5, j+0.5, k+0.5)

    // u-:    dy  dz 
    double u_up = 0.5 * (safe_u(i, j + 1, k) + safe_u(i + 1, j + 1, k));
    double u_down = 0.5 * (safe_u(i, j - 1, k) + safe_u(i + 1, j - 1, k));
    double u_front = 0.5 * (safe_u(i, j, k + 1) + safe_u(i + 1, j, k + 1));
    double u_back = 0.5 * (safe_u(i, j, k - 1) + safe_u(i + 1, j, k - 1));

    // v-:    dx  dz 
    double v_right = 0.5 * (safe_v(i + 1, j, k) + safe_v(i + 1, j + 1, k));
    double v_left = 0.5 * (safe_v(i - 1, j, k) + safe_v(i - 1, j + 1, k));
    double v_front = 0.5 * (safe_v(i, j, k + 1) + safe_v(i, j + 1, k + 1));
    double v_back = 0.5 * (safe_v(i, j, k - 1) + safe_v(i, j + 1, k - 1));

    // w-:    dx  dy 
    double w_right = 0.5 * (safe_w(i + 1, j, k) + safe_w(i + 1, j, k + 1));
    double w_left = 0.5 * (safe_w(i - 1, j, k) + safe_w(i - 1, j, k + 1));
    double w_up = 0.5 * (safe_w(i, j + 1, k) + safe_w(i, j + 1, k + 1));
    double w_down = 0.5 * (safe_w(i, j - 1, k) + safe_w(i, j - 1, k + 1));

    //    curl   
    //  =   u
    double dw_dy = (w_up - w_down) / (2.0 * dy);
    double dv_dz = (v_front - v_back) / (2.0 * dz);
    double wx = dw_dy - dv_dz;

    double du_dz = (u_front - u_back) / (2.0 * dz);
    double dw_dx = (w_right - w_left) / (2.0 * dx);
    double wy = du_dz - dw_dx;

    double dv_dx = (v_right - v_left) / (2.0 * dx);
    double du_dy = (u_up - u_down) / (2.0 * dy);
    double wz = dv_dx - du_dy;

    return Vec3(wx, wy, wz);
}

double MACVectorGrid::vorticityMagnitude(int i, int j, int k) const
{
    Vec3 vort = vorticityVector(i, j, k);
    return vort.norm();  //    Vec3   .norm()
}

double MACVectorGrid::maxVorticity() const
{
    Size3 dim = res();
    double maxVort = 0.0;

#pragma omp parallel //for reduction(max:maxVort)
    for (int k = 0; k < dim[2]; ++k) {
        for (int j = 0; j < dim[1]; ++j) {
            for (int i = 0; i < dim[0]; ++i) {
                double vort = vorticityMagnitude(i, j, k);
                if (vort > maxVort) {
                    maxVort = vort;
                }
            }
        }
    }

    return maxVort;
}

// ============================================================================
//   
// ============================================================================

double MACVectorGrid::computeKineticEnergy() const
{
    Size3 dim = res();
    Vec3 spc = spacing();
    double dV = spc[0] * spc[1] * spc[2];  //   
    double rho = 1.0;  //  (   )

    double kineticEnergy = 0.0;

#pragma omp parallel for reduction(+:kineticEnergy)
    for (int k = 0; k < dim[2]; ++k) {
        for (int j = 0; j < dim[1]; ++j) {
            for (int i = 0; i < dim[0]; ++i) {
                //     
                double uc = 0.5 * (u(i, j, k) + u(i + 1, j, k));
                double vc = 0.5 * (v(i, j, k) + v(i, j + 1, k));
                double wc = 0.5 * (w(i, j, k) + w(i, j, k + 1));

                double speedSquared = uc * uc + vc * vc + wc * wc;
                kineticEnergy += 0.5 * rho * speedSquared * dV;
            }
        }
    }

    return kineticEnergy;
}

double MACVectorGrid::computeEnstrophy() const
{
    Size3 dim = res();
    Vec3 spc = spacing();
    double dV = spc[0] * spc[1] * spc[2];  //   

    double enstrophy = 0.0;

#pragma omp parallel for reduction(+:enstrophy)
    for (int k = 0; k < dim[2]; ++k) {
        for (int j = 0; j < dim[1]; ++j) {
            for (int i = 0; i < dim[0]; ++i) {
                Vec3 vort = vorticityVector(i, j, k);
                double vortMagnitudeSquared = vort[0] * vort[0] +
                    vort[1] * vort[1] +
                    vort[2] * vort[2];
                enstrophy += 0.5 * vortMagnitudeSquared * dV;
            }
        }
    }

    return enstrophy;
}
#endif


