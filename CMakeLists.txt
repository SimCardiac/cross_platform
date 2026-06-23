cmake_minimum_required(VERSION 3.15)
set(PROJECT_NAME Helloworld_PETSc_MPI)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for petsc
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# start a project
project(${PROJECT_NAME} LANGUAGES CXX C)

# Enable CTest for unit testing
enable_testing()

# Find required packages (PETSc should be discoverable via pkg-config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC IMPORTED_TARGET PETSc QUIET)
find_package(Kokkos REQUIRED)
find_package(MPI REQUIRED)

add_library(ins_lib 
    src/poisson/empty.cpp
)
target_include_directories(ins_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ins_lib PUBLIC PkgConfig::PETSC MPI::MPI_CXX)

set(SOLVER_SOURCES
    # Heat equation solvers
    src/heat/heat_DMStag_2D.cpp
    # src/heat/heat_DMStag_CN_2D.cpp
    src/heat/heat_DMStag_edge_2D.cpp
    # Stokes equation solvers
    src/stokes/stokes_DMStag_projection_2D.cpp
    src/stokes/stokes_DMStag_unsteady_2D.cpp
    src/stokes/stokes_DMStag_2D.cpp
    # Unified PDE solver
    src/unified_pde_solver_2D.cpp
    # Nonlinear solvers
    src/nonlinear/test_JFNK.cpp
    src/nonlinear/beam_system.cpp
    # Poisson equation test executables
    src/poisson/tests/poisson_tests_cli.cpp
    src/poisson/tests/poisson_dmstag_tests_cli.cpp
    src/poisson/poisson_DMStag_2D.cpp
    src/poisson/poisson_vertex_2D.cpp
)

foreach(source_file ${SOLVER_SOURCES})

    get_filename_component(exec_name ${source_file} NAME_WE)
    
    add_executable(${exec_name} ${source_file})
    target_include_directories(${exec_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${exec_name} PRIVATE ins_lib)
endforeach()

add_subdirectory(tests/unit/poisson)
add_subdirectory(tests/unit/heat)
add_subdirectory(tests/unit/stokes)

