#include "GridFluidSolver.h"
#include "CenterScalarGrid.h"
#include "Utils.h"

#include <iomanip>
#include <set>
#include <filesystem>
#include <fstream>
#include <string>

#include <vtkStructuredGrid.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkOBJWriter.h>
#include <vtkPolyDataWriter.h>

#define M_PI 3.14159265358979323846

// Funkcije za sejviranje u .vtk formatu
inline double safe_val_vtk(double x);

template<typename CenterScalarGrid, typename MACVectorGrid>
void save_to_vtk(const Size3 Domen, Vec3 res, const std::shared_ptr<CenterScalarGrid> mesh, const std::shared_ptr<CenterScalarGrid> vorticity,
                 const std::shared_ptr<MACVectorGrid> velocity, const std::string& filename);

// Neke globalne promenljive

bool USE_COVECTOR_CORRECTION = true;
bool USE_COLLIDER = false;

int Nx = 80, Ny = 80, Nz = 80;
double dx = 0.05, dy = 0.05, dz = 0.05;

int Scenario = 1;

int numFrames = 100;
int SAVE_VTK = 0;

// Za Scenario 6. (Dva mlaza na susednim stranicama)

double cylinderRadius = 20. * dx;
double cylinderHeight = 20 * dx;
double cylinderOffset = 4. * dx;
double sourceVelocity = 30.;

// MAIN

int main(int argc, char** argv)
{
   // iz komandne linije ...

   if (argc != 12)
   {
       std::cout << "cf3d USE_COVECTOR_CORRECTION USE_COLLIDER Nx Ny Nz dx dy dz Scenario numFrames SAVE_VTK" << std::endl;
       std::cout << "Na primer: cf3d 1 0 80 80 80 0.05 0.05 0.05 1 100 0" << std::endl;
       std::cout << " Scenario 1 : Sfera pada u bazen" << endl;
       std::cout << " Scenario 2 : Taylor-Green vrtlog" << endl;
       std::cout << " Scenario 3 : Pucanje brane" << endl;
       std::cout << " Scenario 4 : Toroidalni torusni vrtlog" << endl;
       std::cout << " Scenario 5 : Poloidalni torusni vrtlog" << endl;
       std::cout << " Scenario 6 : Mlazevi na susednim vertikalnim zidovima" << endl;
       std::cout << " Scenario 7 : Kvadratni vodeni stub u centru" << endl;
       exit(0);
   }

   USE_COVECTOR_CORRECTION = atoi(argv[1]);
   USE_COLLIDER = atoi(argv[2]);
   Nx = atoi(argv[3]);
   Ny = atoi(argv[4]);
   Nz = atoi(argv[5]);
   dx = atof(argv[6]);
   dy = atof(argv[7]);
   dz = atof(argv[8]);
   Scenario = atoi(argv[9]);
   numFrames = atoi(argv[10]);
   SAVE_VTK = atoi(argv[11]);

  // Recompute source/cylinder parameters after parsing command line
  cylinderRadius = (double)Nx / 10.0 * dx;
  cylinderHeight = (double)Nx / 10.0 * dx;

   std::cout << "\nUSE_COVECTOR: " << USE_COVECTOR_CORRECTION << std::endl;
   std::cout << "USE_COLLIDER: " << USE_COLLIDER << std::endl;
   std::cout << "Nx: " << Nx << std::endl;
   std::cout << "Ny: " << Ny << std::endl;
   std::cout << "Nz: " << Nz << std::endl;
   std::cout << "dx: " << dx << std::endl;
   std::cout << "dy: " << dy << std::endl;
   std::cout << "dz: " << dz << std::endl;
   std::cout << "Scenario: " << Scenario << std::endl;
   std::cout << "numFrames: " << numFrames << std::endl;
   std::cout << "SAVE_VTK: " << SAVE_VTK << "\n" << std::endl;

   // Folderi za rezultate

   std::filesystem::path dir_vtk = "results_vtk";
   std::filesystem::path dir_obj = "results_obj";
   if (std::filesystem::exists(dir_vtk)) std::filesystem::remove_all(dir_vtk);
   if (std::filesystem::exists(dir_obj)) std::filesystem::remove_all(dir_obj);

   std::filesystem::create_directory("results_vtk");
   std::filesystem::create_directory("results_obj");

   std::ofstream file("KinetcEnergyAndEnstrophy.txt");

   // Solver

   Size3 Domen(Nx, Ny, Nz);
   Vec3 Origin(0., 0., 0.);
   Vec3 spacing(dx, dy, dz);

   GridFluidSolver solver(Domen, Origin, spacing);

   // Osnovne mreze

   auto free_surface = std::dynamic_pointer_cast<CenterScalarGrid>(solver.grids().getScalarGrid("FREE_SURFACE"));
   auto vorticity_magnitude = std::dynamic_pointer_cast<CenterScalarGrid>(solver.grids().getScalarGrid("VORTICITY"));
   auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solver.grids().getVelocityGrid());

   auto sourceMask = std::dynamic_pointer_cast<CenterScalarGrid>(solver.grids().getScalarGrid("SOURCE_MASK"));
   sourceMask->fillData(1.0);  // inicijalizuj na 1 (nije source)

   if(Scenario == 1) 
   {
    // Scenario 1 - Sfera pada u bazen

    const double radius = 0.25 * Domen[0] * spacing[0];
    Vec3 sphereCenter(0.5 * spacing[0] * Domen[0], 0.7 * spacing[1] * Domen[1], 0.5 * spacing[2] * Domen[2]);

    free_surface->fillData([&](double x, double y, double z) -> double
        {
            // Izračunaj kolajder SDF (ako postoji)
            if(USE_COLLIDER)
            {
             double colliderSDF = solver.colliderSDF->sample(Vec3(x,y,z));
             if(colliderSDF < 0)
             {
				 return 2; // unutar kolajdera - nema fluida
             }
            }

            // 1. SDF za sferu: rastojanje od tacke do povrsine sfere
            double dx = x - sphereCenter[0];
            double dy = y - sphereCenter[1];
            double dz = z - sphereCenter[2];
            double distToCenter = std::sqrt(dx * dx + dy * dy + dz * dz);
            double sdfSphere = distToCenter - radius;

            // 2. SDF za polubeskonacnu vodu (poluprostor y <= waterSurfaceY)
            double sdfWater = y - 0.15 * spacing[1] * Domen[1];

            // 3. UNION: najblizi interfejs odredjuje SDF vrednost
            return (sdfSphere < sdfWater) ? sdfSphere : sdfWater;
        });


   }

   else if(Scenario == 2)
   {
    // Scenario 2 - Taylor-Green vrtlog

    //Vec3 spacing(0.157079632, 0.157079632, 0.157079632); ili 0.314159264
    //Size3 Domen(80, 100, 80);

    //double Lx = 2 * M_PI, Ly = 2 * M_PI, Lz = 2 * M_PI;
    double dx = spacing[0], dy = spacing[1], dz = spacing[2];

    free_surface->fillData([&](double x, double y, double z) -> double
        {
            if(y > 0.8 * Domen[1] * dy) return 5;
            else return -5;
        });

    double G = 50.;

    velocity->forEachU([&](int i, int j, int k) {
             if (j * dy < 0.8 * Domen[1] * dy)
             {
              velocity->setUDataAt(i, j, k, G * std::sin(i * dx) * std::cos(j * dy) * std::cos(k * dz));
              //velocity->setUDataAt(i + 1, j, k, G * std::sin((i + 1) * dx) * std::cos(j * dy) * std::cos(k * dz));
             }
    });

    velocity->forEachW([&](int i, int j, int k) {
        if (j * dy < 0.8 * Domen[1] * dy)
        {
            velocity->setWDataAt(i, j, k, -G * std::cos(i * dx) * std::cos(j * dy) * std::sin(k * dz));
            //velocity->setWDataAt(i, j, k + 1, -G * std::cos(i * dx) * std::cos(j * dy) * std::sin((k + 1) * dz));
        }
    });

   }

   else if(Scenario == 3)
   {
    // Scenario 3 - Pucanje brane

    double waterWidth = 0.3 * Domen[0] * spacing[0];   // voda zauzima 30% sirine u X
    double waterHeight = 0.7 * Domen[1] * spacing[1];   // voda zauzima 70% visine u Y

    free_surface->fillData([&](double x, double y, double z) -> double
        {
            // SDF za vodu: negativno unutar vode, pozitivno spolja
            // Voda je u regionu: x < waterWidth AND y < waterHeight
            double sdf_x = x - waterWidth; 
            double sdf_y = y - waterHeight;

            // Intersection: voda je gde su OBA uslova zadovoljena
            double sdfWater = std::max(sdf_x, sdf_y);

            if(y <= 0.1 * waterHeight && x >= waterWidth) return -y;
            else return sdfWater;
        });
    solver.reinitializeFreeSurface();
   }

   else if(Scenario == 4)
   {
    // Scenario 4 - Toroidalni vrtlozni prsten

    double water_level = 0.6;                             // relativni nivo vode
    double cx = dx * Domen[0] / 2., 
           cy = water_level * dy * Domen[1] / 2, 
           cz = dz * Domen[2] / 2.;                       // centar torusa

    double R = 0.3 * Domen[0] * dx;                       // veliki radijus (rastojanje od centra do obruca)
    double r = 0.1 * Domen[0] * dx;                       // mali radijus (poprecni presek)
    double circulation = 10.;                             // jacina vrtloga

    std::cout << "Toroidalni vrtlog:  R = " << R << "   r = " << r << std::endl;

    free_surface->fillData([&](double x, double y, double z) -> double
        {
            if (y > water_level * Domen[1] * dy) return 5;
            else return -5;
        });

    velocity->forEachU([&](int i, int j, int k) {
        if (j * dy < water_level * Domen[1] * dy)
        {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dy;
            double z = (k + 0.5) * dz;

            // Vektor od centra torusa u (cx, cy, cz) do trenutne tacke, bez y komponente
            double px = x - cx;
            double pz = z - cz;
            double radial_dist = sqrt(px * px + pz * pz);  // rastojanje u XZ ravni

            // Projekcija tacke na centralni krug u ravni normalnoj na Y
            double proj_x, proj_z;
            if (radial_dist > 1e-8) {
                proj_x = cx + R * px / radial_dist;
                proj_z = cz + R * pz / radial_dist;
            }
            else {
                proj_x = cx + R;
                proj_z = cz;
            }
            double proj_y = cy;

            // Vektor od projekcije do tacke (udaljenost od prstena)
            double vx = x - proj_x;
            double vy = y - cy;
            double vz = z - proj_z;
            double dist_to_ring = sqrt(vx * vx + vy * vy + vz * vz);

            // Tangenta na torus: rotacija za 90 stepeni u XZ ravni (tangenta na krug)
            // Za Y-osu: tangenta je ( -dz, 0, dx ) normalizovano

            double tangent_x = -pz;
            double tangent_y = 0.0;
            double tangent_z = px;

            double len = sqrt(tangent_x * tangent_x + tangent_z * tangent_z);
            if (len > 1e-8) {
                tangent_x /= len;
                tangent_z /= len;
            }

            // Brzina unutar poprecnog preseka
            if (dist_to_ring < r) {
                double core_ratio = (r - dist_to_ring) / r;
                if (core_ratio < 0) core_ratio = 0;
                double speed = (circulation / (2.0 * M_PI * R)) * core_ratio;
                velocity->setUDataAt(i, j, k, speed * tangent_x);
            }
            else {
                // Spolja: brzina opada do nule
                double decay = std::max(0.0, (r + 0.5 * dx - dist_to_ring) / (0.5 * dx));
                double speed = (circulation / (2.0 * M_PI * R)) * decay;
                velocity->setUDataAt(i, j, k, speed * tangent_x);
            }
        }});

    velocity->forEachV([&](int i, int j, int k) {
        if (j * dy < water_level * Ny * dy)
        {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dy;
            double z = (k + 0.5) * dz;

            // Vektor od centra torusa u (cx, cy, cz) do trenutne take, bez y komponente
            double px = x - cx;
            double pz = z - cz;
            double radial_dist = sqrt(px * px + pz * pz);  // rastojanje u XZ ravni

            // Projekcija take na centralni krug u ravni normalnoj na Y
            double proj_x, proj_z;
            if (radial_dist > 1e-8) {
                proj_x = cx + R * px / radial_dist;
                proj_z = cz + R * pz / radial_dist;
            }
            else {
                proj_x = cx + R;
                proj_z = cz;
            }
            double proj_y = cy;

            // Vektor od projekcije do take (udaljenost od prstena)
            double vx = x - proj_x;
            double vy = y - cy;
            double vz = z - proj_z;
            double dist_to_ring = sqrt(vx * vx + vy * vy + vz * vz);

            // Tangenta na torus: rotacija za 90 u XZ ravni (tangenta na krug)
            // Za Y-osu: tangenta je ( -dz, 0, dx ) normalizovano

            double tangent_x = -pz;
            double tangent_y = 0.0;
            double tangent_z = px;

            double len = sqrt(tangent_x * tangent_x + tangent_z * tangent_z);
            if (len > 1e-8) {
                tangent_x /= len;
                tangent_z /= len;
            }

            // Brzina unutar poprenog preseka
            if (dist_to_ring < r) {
                double core_ratio = (r - dist_to_ring) / r;
                if (core_ratio < 0) core_ratio = 0;
                double speed = (circulation / (2.0 * M_PI * R)) * core_ratio;
                velocity->setVDataAt(i, j, k, speed * tangent_y);
            }
            else {
                // Spolja: brzina opada do nule
                double decay = std::max(0.0, (r + 0.5 * dx - dist_to_ring) / (0.5 * dx));
                double speed = (circulation / (2.0 * M_PI * R)) * decay;
                velocity->setVDataAt(i, j, k, speed * tangent_y);
            }
        }});

    velocity->forEachW([&](int i, int j, int k) {
        if (j * dy < water_level * Ny * dy)
        {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dy;
            double z = (k + 0.5) * dz;

            // Vektor od centra torusa u (cx, cy, cz) do trenutne take, bez y komponente
            double px = x - cx;
            double pz = z - cz;
            double radial_dist = sqrt(px * px + pz * pz);  // rastojanje u XZ ravni

            // Projekcija take na centralni krug u ravni normalnoj na Y
            double proj_x, proj_z;
            if (radial_dist > 1e-8) {
                proj_x = cx + R * px / radial_dist;
                proj_z = cz + R * pz / radial_dist;
            }
            else {
                proj_x = cx + R;
                proj_z = cz;
            }
            double proj_y = cy;

            // Vektor od projekcije do take (udaljenost od prstena)
            double vx = x - proj_x;
            double vy = y - cy;
            double vz = z - proj_z;
            double dist_to_ring = sqrt(vx * vx + vy * vy + vz * vz);

            // Tangenta na torus: rotacija za 90 u XZ ravni (tangenta na krug)
            // Za Y-osu: tangenta je ( -dz, 0, dx ) normalizovano

            double tangent_x = -pz;
            double tangent_y = 0.0;
            double tangent_z = px;

            double len = sqrt(tangent_x * tangent_x + tangent_z * tangent_z);
            if (len > 1e-8) {
                tangent_x /= len;
                tangent_z /= len;
            }

            // Brzina unutar poprenog preseka
            if (dist_to_ring < r) {
                double core_ratio = (r - dist_to_ring) / r;
                if (core_ratio < 0) core_ratio = 0;
                double speed = (circulation / (2.0 * M_PI * R)) * core_ratio;
                velocity->setWDataAt(i, j, k, speed * tangent_z);
            }
            else {
                // Spolja: brzina opada do nule
                double decay = std::max(0.0, (r + 0.5 * dx - dist_to_ring) / (0.5 * dx));
                double speed = (circulation / (2.0 * M_PI * R)) * decay;
                velocity->setWDataAt(i, j, k, speed * tangent_z);
            }
        }});
    }

else if (Scenario == 5)
{
    // Scenario 5 - Polidalni vrtlozni prsten

    double water_level = 0.6;
    double cx = dx * Domen[0] / 2.,
        cy = water_level * dy * Domen[1] / 2,
        cz = dz * Domen[2] / 2.;

    double R = 0.3 * Domen[0] * dx;   // Veliki poluprecnik
    double a = 0.1 * Domen[0] * dx;   // Poluprecnik jezgra
    double circulation = 10.;         // Gamma
    int    direction = -1;          // +1 ili -1 za kontrolu smera (+Y ili -Y)

    std::cout << "Vortex Ring: R = " << R << "  a = " << a << "  Gamma = " << circulation << std::endl;

    // Slobodna povrsina
    free_surface->fillData([&](double x, double y, double z) -> double {
        return (y > water_level * Domen[1] * dy) ? 5.0 : -5.0;
        });

    // Lambda za brzinu - Racuna za SVAKU tacku u fluidu
    auto getVortexVelocity = [&](double x, double y, double z, double& speed_out) -> std::array<double, 3> {
        double px = x - cx, pz = z - cz;
        double rho = sqrt(px * px + pz * pz);

        // Na samoj Y-osi: brzina je nula (singularitet koordinata)
        if (rho < 1e-8) {
            speed_out = 0.0;
            return { 0.0, 0.0, 0.0 };
        }

        // Koordinate u meridionalnoj (rho, y) ravni
        double drho = rho - R;
        double dy = y - cy;
        double d2 = drho * drho + dy * dy;
        double d = sqrt(d2);

        // Rankine profil: cvrsta rotacija unutra, potencijalno polje spolja
        double v_theta;
        if (d <= a) {
            // Unutar jezgra: linearni rast od 0 u centru do max na ivici
            v_theta = (circulation * d) / (2.0 * M_PI * a * a);
        }
        else {
            // Spolja: opadanje ~1/d za dugodometnu indukciju
            v_theta = circulation / (2.0 * M_PI * d);
        }

        // Tangenta za CCW cirkulaciju u (rho,y) ravni x smer
        double t_rho = -dy / (d + 1e-8) * direction;
        double t_y = drho / (d + 1e-8) * direction;

        // Mapiranje u 3D Dekartove komponente
        double u = v_theta * t_rho * (px / rho);
        double v = v_theta * t_y;
        double w = v_theta * t_rho * (pz / rho);

        // Magnituda za dijagnostiku
        speed_out = sqrt(u * u + v * v + w * w);
        return { u, v, w };
        };

    // Inicijalizacija na MAC stranicama - BEZ uslova koji iskljucuju unutrasnjost!
    velocity->forEachU([&](int i, int j, int k) {
        // Provera da li smo u fluidu (ispod slobodne povrine)
        if (j * dy < water_level * Domen[1] * dy) {
            double x = (i + 0.5) * dx, y = (j + 0.5) * dy, z = (k + 0.5) * dz;
            double spd;
            auto vel = getVortexVelocity(x, y, z, spd);
            velocity->setUDataAt(i, j, k, vel[0]);
        }
        });
    velocity->forEachV([&](int i, int j, int k) {
        if (j * dy < water_level * Domen[1] * dy) {
            double x = (i + 0.5) * dx, y = (j + 0.5) * dy, z = (k + 0.5) * dz;
            double spd;
            auto vel = getVortexVelocity(x, y, z, spd);
            velocity->setVDataAt(i, j, k, vel[1]);
        }
        });
    velocity->forEachW([&](int i, int j, int k) {
        if (j * dy < water_level * Domen[1] * dy) {
            double x = (i + 0.5) * dx, y = (j + 0.5) * dy, z = (k + 0.5) * dz;
            double spd;
            auto vel = getVortexVelocity(x, y, z, spd);
            velocity->setWDataAt(i, j, k, vel[2]);
        }
        });

    // Dijagnostika: Ispis statistike za ceo presek torusa
    double min_spd = 1e9, max_spd = 0.0, avg_spd = 0.0;
    int count_core = 0, count_total = 0;

    for (int i = 0; i < Domen[0]; ++i)
        for (int j = 0; j < Domen[1]; ++j)
            for (int k = 0; k < Domen[2]; ++k) {
                if (j * dy >= water_level * Domen[1] * dy) continue; // van fluida

                double x = (i + 0.5) * dx, y = (j + 0.5) * dy, z = (k + 0.5) * dz;
                double px = x - cx, pz = z - cz;
                double rho = sqrt(px * px + pz * pz);
                if (rho < 1e-8) continue;

                double drho = rho - R, dy = y - cy;
                double d = sqrt(drho * drho + dy * dy);

                if (d <= a + 2.0 * dx) { // zona uticaja (jezgro + tranzicija)
                    count_total++;
                    double u = velocity->u(i, j, k);
                    double v = velocity->v(i, j, k);
                    double w = velocity->w(i, j, k);
                    double spd = sqrt(u * u + v * v + w * w);

                    if (d <= a) { // striktno unutar jezgra
                        count_core++;
                        avg_spd += spd;
                        if (spd < min_spd) min_spd = spd;
                        if (spd > max_spd) max_spd = spd;
                    }
                }
            }

    if (count_core > 0) avg_spd /= count_core;
    std::cout << "\n=== Dijagnostika torusa ===" << std::endl;
    std::cout << "Celije u jezgru (d<a): " << count_core << std::endl;
    std::cout << "Celije u zoni uticaja: " << count_total << std::endl;
    std::cout << "Brzina u jezgru: min=" << min_spd << "  avg=" << avg_spd
        << "  max=" << max_spd << std::endl;
    std::cout << "Teorijski max (Gamma/2 Pi a): " << (circulation / (2 * M_PI * a)) << std::endl;
    std::cout << "========================\n" << std::endl;
    }

else if (Scenario == 6)
{
    // === Scenario 6: Two cylindrical sources on adjacent walls ===

    //cylinderRadius = 12. * dx;
    //cylinderHeight = 10. * dx;
    //cylinderOffset = 4. * dx;

    double Lx = Domen[0] * dx;
    double Ly = Domen[1] * dy;
    double Lz = Domen[2] * dz;

    // === SOURCE 1: Left wall (X=0), cylinder offset into domain ===
    // Base sits at x = cylinderOffset (not x=0)
    // Center is half the cylinder height from the base
    Vec3 cyl1_base(cylinderOffset, 0.8 * Ly, 0.7 * Lz);
    Vec3 cyl1_center(cylinderOffset + cylinderHeight / 2.0, 0.8 * Ly, 0.7 * Lz);

    // === SOURCE 2: Back wall (Z=Lz), cylinder offset into domain ===
    // Base sits at z = Lz - cylinderOffset
    Vec3 cyl2_base(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset);
    Vec3 cyl2_center(0.7 * Lx, 0.8 * Ly, Lz - cylinderOffset - cylinderHeight / 2.0);

    // === INITIALIZE SDF ===
    free_surface->fillData([&](double x, double y, double z) -> double
        {
            // --- CYLINDER 1: Axis along +X, base at x = cylinderOffset ---
            double dy1 = y - cyl1_center[1];
            double dz1 = z - cyl1_center[2];
            double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);

            // Aksijalno rastojanje: mereno od baze (cylinderOffset)
            double axial_dist1 = x - cylinderOffset;

            bool inside1 = (radial_dist1 < cylinderRadius &&
                axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight);

            // --- CYLINDER 2: Axis along -Z, base at z = Lz - cylinderOffset ---
            double dx2 = x - cyl2_center[0];
            double dy2 = y - cyl2_center[1];
            double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);

            // Aksijalno rastojanje: mereno od baze (Lz - cylinderOffset)
            double axial_dist2 = (Lz - cylinderOffset) - z;

            bool inside2 = (radial_dist2 < cylinderRadius &&
                axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight);

            // Explicit: negative inside (fluid), positive outside (air)
            if (inside1 || inside2) {
                return -1.0;  // fluid (strong negative for stability)
            }
            else {
                return 5.0;   // air
            }
        });

    // solver.reinitializeFreeSurface();

    // === Dijagnostika ===

    std::cout << "\n=== Provera SDF ===" << std::endl;
    int test_i = (int)((cylinderOffset + cylinderHeight / 2.0) / dx);
    int test_j = (int)(0.8 * Ly / dy);  // Y centar cilindra 1
    int test_k = (int)(0.7 * Lz / dz);  // Z centar cilindra 1

    if (test_i >= 0 && test_i < Domen[0] &&
        test_j >= 0 && test_j < Domen[1] &&
        test_k >= 0 && test_k < Domen[2]) {

        double sdf_val = free_surface->at(test_i, test_j, test_k);
        std::cout << "Test tacka: i=" << test_i << ", j=" << test_j << ", k=" << test_k << std::endl;
        std::cout << "Fizicki polozaj: x=" << (test_i + 0.5) * dx
            << ", y=" << test_j * dy
            << ", z=" << test_k * dz << std::endl;
        std::cout << "SDF u centru cilindra 1: " << sdf_val << std::endl;
        std::cout << "cylinderOffset = " << cylinderOffset << std::endl;
        std::cout << "cylinderHeight = " << cylinderHeight << std::endl;
        if (sdf_val < 0.0) {
            std::cout << "OK: SDF je negativan - fluid!" << std::endl;
        }
        else {
            std::cout << "ERROR: SDF je pozitivan - nije fluid!" << std::endl;
        }
    }
    else {
        std::cout << "ERROR: Test indeksi su izvan opsega mreze!" << std::endl;
    }
    std::cout << "====================\n" << std::endl;

    // === INITIALIZE VELOCITY ===
    velocity->fillData([&](double x, double y, double z) -> Vec3
        {
            // Cylinder 1
            double dy1 = y - cyl1_center[1];
            double dz1 = z - cyl1_center[2];
            double radial_dist1 = std::sqrt(dy1 * dy1 + dz1 * dz1);
            double axial_dist1 = x - cylinderOffset;  // measured from offset base

            if (radial_dist1 < cylinderRadius &&
                axial_dist1 >= 0.0 && axial_dist1 <= cylinderHeight) {
                return Vec3(sourceVelocity, 0.0, 0.0);  // +X direction
            }

            // Cilindar 2
            double dx2 = x - cyl2_center[0];
            double dy2 = y - cyl2_center[1];
            double radial_dist2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
            double axial_dist2 = (Lz - cylinderOffset) - z;  //     

            if (radial_dist2 < cylinderRadius &&
                axial_dist2 >= 0.0 && axial_dist2 <= cylinderHeight) {
                return Vec3(0.0, 0.0, -sourceVelocity);  // -Z 
            }    

            return Vec3(0.0, 0.0, 0.0);
        });

        std::cout << "\n=== Scenario 6: Two cylindrical sources (with offset) ===" << std::endl;
        std::cout << "  Cylinder 1: Base at X=" << cylinderOffset << ", extends to X=" << (cylinderOffset + cylinderHeight) << std::endl;
        std::cout << "  Cylinder 2: Base at Z=" << (Lz - cylinderOffset) << ", extends to Z=" << (Lz - cylinderOffset - cylinderHeight) << std::endl;
        std::cout << "  Radius: " << cylinderRadius << "m, Height: " << cylinderHeight << "m" << std::endl;
        std::cout << "  Offset from wall: " << cylinderOffset << "m (" << (cylinderOffset / dx) << " cells)" << std::endl;
    }

    else if(Scenario == 7)
    {
        // Scenario 7 - Kvadratni vodeni stub u centru

        double waterHeight = 0.9 * Domen[1] * dy;   // voda zauzima 90% visine u Y
        double cx = 0.5 * Domen[0] * dx;
        double cz = 0.5 * Domen[2] * dz;
        double halfWidth = 0.1 * std::min(Domen[0] * dx, Domen[2] * dz);

        free_surface->fillData([&](double x, double y, double z) -> double {
            // Kvadrat u X-Z ravni
            double sdf_square = std::max(std::abs(x - cx) - halfWidth,
                                std::abs(z - cz) - halfWidth);
            double sdf_water = y - waterHeight;

            if(y <= 0.1 * waterHeight) return -2;
            else return std::max(sdf_square, sdf_water);
            });
        solver.reinitializeFreeSurface();
    }

    // ====================== Simulation start=========================

    if (SAVE_VTK == 1)
    {
          vorticity_magnitude->fillData([&](int i, int j, int k) -> double
            {
                return 0.;  //velocity->vorticityMagnitude(i, j, k);
            });
    }

    if (SAVE_VTK == 1) save_to_vtk(Domen, spacing, free_surface, vorticity_magnitude, velocity, "results_vtk/fluid_0.vtk");

    auto ukupno_vreme_start =  std::chrono::high_resolution_clock::now();

    for (int frame = 1; frame <= numFrames; frame++)
    {
        auto startTime = std::chrono::high_resolution_clock::now();

        //if(frame == 1 && Scenario == 6) solver.resetSourceVelocity(cylinderRadius, cylinderHeight, cylinderOffset, sourceVelocity);

        // Update values
        std::cout << "calculating frame " << frame << "\n" << std::endl;
        solver.update(frame - 1);
        std::cout.clear();

        // Totalna kineticka enegija i enstrofija

        //double kinetic_energy = velocity->computeKineticEnergy();
        //double enstrophy = velocity->computeEnstrophy();

        //std::cout << "KE = " << kinetic_energy << " Enstrophy = " << enstrophy << "\n" << std::endl;
        //file << "KE = " << kinetic_energy << " Enstrophy = " << enstrophy << std::endl;

        // Save results

        std::cout << "Saving results ...\n\n" << std::endl;
        if(frame % 5 == 0 || frame == 1 || frame == numFrames)
        {
         if(SAVE_VTK == 1)
         {
             // Popunite mrezu vorticiteta koristeci novu metodu
             vorticity_magnitude->fillData([&](int i, int j, int k) -> double
                 {
                     return 0.;  //velocity->vorticityMagnitude(i, j, k);
                 });
         }

         // Izdvajanje nultog nivoa - slobodne povrsine

         vtkSmartPointer<vtkStructuredGrid> ResultsDataSet = vtkSmartPointer<vtkStructuredGrid>::New();
         ResultsDataSet->SetDimensions(Domen[0], Domen[1], Domen[2]);

         // Kreiranje tacka
         vtkNew<vtkPoints> points;
         points->SetNumberOfPoints(Domen[0] * Domen[1] * Domen[2]);

         // Kreiranje niza za podatke
         vtkNew<vtkFloatArray> dataArray;
         dataArray->SetName("vrednosti");
         dataArray->SetNumberOfValues(Domen[0] * Domen[1] * Domen[2]);

         int pointId = 0;
         for (int k = 0; k < Domen[2]; k++) {
            for (int j = 0; j < Domen[1]; j++) {
                for (int i = 0; i < Domen[0]; i++) {
                    
                    double x = (i - 1) * spacing[0];
                    double y = (j - 1) * spacing[1];
                    double z = (k - 1) * spacing[2];

                    points->SetPoint(pointId, x, y, z);
                    dataArray->SetValue(pointId, free_surface->at(i, j, k));
                    pointId++;
                }
            }
         }

         ResultsDataSet->SetPoints(points);
         ResultsDataSet->GetPointData()->SetScalars(dataArray);

         vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
         contourFilter->SetInputData(ResultsDataSet);
         contourFilter->GenerateValues(1, 0, 0);
         contourFilter->Update();

         // zapis u .obj
         vtkNew<vtkOBJWriter> writer;
         std::stringstream file_name;
         file_name << "results_obj/frame_" << std::setfill('0') << std::setw(6) << frame << ".obj";
         writer->SetFileName(file_name.str().c_str());
         writer->SetInputConnection(contourFilter->GetOutputPort());
         writer->Write();

         if(SAVE_VTK == 1) save_to_vtk(Domen, spacing, free_surface, vorticity_magnitude, velocity, "results_vtk/fluid_" + std::to_string(frame) + ".vtk");

         auto stopTime = std::chrono::high_resolution_clock::now();
         double dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
         std::cout << "Frame: " << frame << "  finished in  " << dt << " sec\n" << std::endl;

        }
    }
    file.close();

    auto ukupno_vreme_end =  std::chrono::high_resolution_clock::now();
	double ukupno_vreme = std::chrono::duration<float, std::chrono::seconds::period>(ukupno_vreme_end - ukupno_vreme_start).count();
	std::cout << "Ukupno vreme simulacije: " << ukupno_vreme << " sec." << endl;

    return 0;
}

// Zapis u .vtk formatu

inline double safe_val_vtk(double x) {
    if (std::isnan(x) || std::isinf(x)) return 0.0;
    if (x > 1e6) return 1e6;
    if (x < -1e6) return -1e6;
    if (fabs(x) < 1e-6) return 0;
    return x;
}

template<typename CenterScalarGrid, typename MACVectorGrid>
void save_to_vtk(const Size3 Domen, Vec3 res, const std::shared_ptr<CenterScalarGrid> mesh, const std::shared_ptr<CenterScalarGrid> vorticity,
    const std::shared_ptr<MACVectorGrid> velocity, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Ne mogu da otvorim fajl za upis: " + filename);
    }

    const int nx = Domen[0];
    const int ny = Domen[1];
    const int nz = Domen[2];

    // Podrazumevane vrednosti za mreu
    double origin[3] = { 0.0, 0.0, 0.0 };
    double spacing[3] = { res[0], res[1], res[2] };  // dx, dy, dz

    // --- VTK zaglavlje ---
    file << "# vtk DataFile Version 4.0\n";
    file << "3D Scalar Field from C++ Grid\n";
    file << "ASCII\n";  // Moe i BINARY, ali ASCII je laki za debagovanje
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "ORIGIN " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    file << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << "\n";
    file << "POINT_DATA " << nx * ny * nz << "\n";

    file << "SCALARS " << "surface" << " double 1\n";
    file << "LOOKUP_TABLE default\n";

    // slobodna povrsina

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                file << safe_val_vtk(mesh->at(i, j, k)) << "\n";
            }
        }
    }

    file << "SCALARS " << "vorticity" << " double 1\n";
    file << "LOOKUP_TABLE default\n";

    // vrtloznost

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                file << safe_val_vtk(vorticity->at(i, j, k)) << "\n";
            }
        }
    }

    // brzine u, v i w

    file << "VECTORS velocity double\n";
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // Interpolacija u centru elije
                double ui = 0.5 * (velocity->u(i, j, k) + velocity->u(i + 1, j, k));
                double vj = 0.5 * (velocity->v(i, j, k) + velocity->v(i, j + 1, k));
                double wk = 0.5 * (velocity->w(i, j, k) + velocity->w(i, j, k + 1));
                file << safe_val_vtk(ui) << " " << safe_val_vtk(vj) << " " << safe_val_vtk(wk) << "\n";
            }
        }
    }

    file.close();
}


