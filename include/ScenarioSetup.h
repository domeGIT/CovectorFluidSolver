#pragma once

#include "GridFluidSolver.h"
#include "CenterScalarGrid.h"
#include "MACVectorGrid.h"

// Setup scenario-specific initial conditions. This function mirrors the
// scenario blocks previously located inside main.cpp.
void setupScenario(const Size3& Domen, const Vec3& spacing, GridFluidSolver& solver,
				   std::shared_ptr<CenterScalarGrid> free_surface,
				   std::shared_ptr<CenterScalarGrid> vorticity_magnitude,
				   std::shared_ptr<MACVectorGrid> velocity);
