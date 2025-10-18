cmake_minimum_required(VERSION 3.15)
set(PROJECT_NAME Helloworld_PETSc_MPI)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for petsc
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

# start a project
project(${PROJECT_NAME} LANGUAGES CXX C)
add_executable(ex1 src/main.cpp)
add_executable(ex2 src/DMStag_2D.cpp)
add_executable(ex3 
    src/wyh/oscillating_disk/wyh_solve_fluid_structure_interaction_system3D200630.cpp
    src/wyh/oscillating_disk/wyh_quadratic_rectangular_element_interpolation3D200527.cpp
    src/wyh/oscillating_disk/main_fluid_solid200508.cpp
    src/wyh/oscillating_disk/wyh_output_result.cpp
)

# find packages
include(FindPkgConfig REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Kokkos REQUIRED)
find_package(MPI REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

# link libraries
target_link_libraries(ex1 PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
target_link_libraries(ex2 PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
target_link_libraries(ex3 PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)

target_include_directories(ex1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# io: parallel VTI writer
add_library(io_vti src/io/vti.cpp)
target_include_directories(io_vti PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io)
target_link_libraries(io_vti PUBLIC MPI::MPI_CXX)

# Example: parallel VTI output
add_executable(ex_vti src/examples/vti_example.cpp)
target_link_libraries(ex_vti PRIVATE io_vti MPI::MPI_CXX)

# grid: space-filling curve ordering library
add_library(grid_sfc src/grid/space_filling_order.cpp)
target_include_directories(grid_sfc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/grid)

# Example: space-filling curve test
add_executable(ex_sfc src/examples/sfc_test.cpp)
target_link_libraries(ex_sfc PRIVATE grid_sfc)

# Example: space-filling curve with PETSc DMDA
add_executable(ex_sfc_dmda src/examples/sfc_dmda_example.cpp)
target_link_libraries(ex_sfc_dmda PRIVATE grid_sfc PkgConfig::PETSC MPI::MPI_CXX)
