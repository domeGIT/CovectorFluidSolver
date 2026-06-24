cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
  cmake_policy(SET CMP0071 NEW)
endif()

project(CovectorFluid3D)

SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_STANDARD 17)

# Target name
set(target cf3d)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Pronađi OpenMP
find_package(OpenMP REQUIRED)

# Includes
include_directories(
.
include
D:/Racun/eigen-3.4.0/Eigen3/include)

LINK_DIRECTORIES()

# Sources
set(sources
		src/ScalarGrid.cpp 
		src/CollocatedGrid.cpp
        src/MACVectorGrid.cpp 
		src/Utils.cpp
        src/ScalarBuilders.cpp 
		src/VectorBuilders.cpp
        src/GridManager.cpp 
        src/GridFluidSolver.cpp 
		src/FluidSolver.cpp 
		src/PressureSolver.cpp
		src/SemiLagrangianAdvectSolver.cpp 
		src/main.cpp) 
		
set(headers
        include/CenterScalarGrid.h 
		include/VertextScalarGrid.h 
		include/CollocatedGrid.h 
		include/CenterVectorGrid.h 
		include/VertexVectorGrid.h
        include/MACVectorGrid.h 
		include/Utils.h 
        include/BaseGrid.h 
		include/ScalarGrid.h 
		include/VectorGrid.h
		include/GridManager.h 
		include/ScalarBuilders.h
        include/VectorBuilders.h 
		include/GridFluidSolver.h
		include/FluidSolver.h 
        include/PressureSolver.h 
		include/SemiLagrangianAdvectSolver.h)

# Build executable
add_executable(${target} ${sources}  ${headers} src/main.cpp)
 
set(EXECUTABLE_OUTPUT_PATH ../bin/)

# Ensure target has absolute include directory to avoid relative include issues
target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(${target}
        PRIVATE
        ${DEFAULT_LINKER_OPTIONS}
		${VTK_LIBRARIES}
		OpenMP::OpenMP_CXX)
