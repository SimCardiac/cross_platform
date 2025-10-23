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

# Enable CTest for unit testing
enable_testing()
add_executable(ex1 src/main.cpp)
add_executable(ex_steady_stokes_stagger src/steady_stokes_DMStag_2D.cpp)
add_executable(ex_poisson_stagger src/poisson_DMStag_2D.cpp)
add_executable(ex_poisson_center src/poisson_DM_2D.cpp)
add_executable(ex_heat_stagger src/heat_DMStag_2D.cpp src/DMStag_boundary_helpers.cpp)
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
target_link_libraries(ex_poisson_center PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
target_link_libraries(ex_poisson_stagger PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
target_link_libraries(ex_heat_stagger PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
target_link_libraries(ex_steady_stokes_stagger PkgConfig::PETSC MPI::MPI_CXX)
target_link_libraries(ex3 PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)

target_include_directories(ex1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex_poisson_center PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex_poisson_stagger PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex_heat_stagger PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex_steady_stokes_stagger PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(ex3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# io: parallel VTI writer
add_library(io_vti src/io/vti.cpp)
target_include_directories(io_vti PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io)
target_link_libraries(io_vti PUBLIC MPI::MPI_CXX)

# Example: parallel VTI output
add_executable(ex_vti src/examples/vti_example.cpp)
target_link_libraries(ex_vti PRIVATE io_vti MPI::MPI_CXX)

# Poisson testing infrastructure
add_library(poisson_mms src/poisson/poisson_mms.cpp)
target_include_directories(poisson_mms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(poisson_metrics src/poisson/poisson_metrics.cpp)
target_link_libraries(poisson_metrics PUBLIC PkgConfig::PETSC MPI::MPI_CXX)
target_include_directories(poisson_metrics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(poisson_solver src/poisson/poisson_solver.cpp)
target_link_libraries(poisson_solver PUBLIC PkgConfig::PETSC MPI::MPI_CXX)
target_include_directories(poisson_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Poisson CLI test runner
add_executable(poisson_tests_cli src/poisson/tests/poisson_tests_cli.cpp)
target_link_libraries(poisson_tests_cli PRIVATE 
    poisson_mms poisson_metrics poisson_solver 
    PkgConfig::PETSC MPI::MPI_CXX)
target_include_directories(poisson_tests_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add tests subdirectory
add_subdirectory(tests/unit/poisson)

