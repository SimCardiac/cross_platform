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
