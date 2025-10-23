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

# # Poisson libraries (DMDA-based and DMStag-based)
# foreach(lib_name poisson_mms poisson_metrics poisson_solver poisson_dmstag_metrics poisson_dmstag_solver)
#     add_library(${lib_name} src/poisson/${lib_name}.cpp)
#     target_include_directories(${lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
#     target_link_libraries(${lib_name} PUBLIC PkgConfig::PETSC MPI::MPI_CXX)
# endforeach()

# INS library combining all poisson modules
add_library(ins_lib 
    src/poisson/poisson_mms.cpp
    src/poisson/poisson_metrics.cpp
    src/poisson/poisson_solver.cpp
    src/poisson/poisson_dmstag_metrics.cpp
    src/poisson/poisson_dmstag_solver.cpp
)
target_include_directories(ins_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ins_lib PUBLIC PkgConfig::PETSC MPI::MPI_CXX)

# Poisson executables (CLI test runners and examples)
foreach(exec_name poisson_tests_cli poisson_dmstag_tests_cli poisson_DMStag_2D)
    add_executable(${exec_name} src/poisson/tests/${exec_name}.cpp)
    target_include_directories(${exec_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${exec_name} PRIVATE PkgConfig::PETSC MPI::MPI_CXX ins_lib)
endforeach()

# Fix path for poisson_DMStag_2D (not in tests/ subdirectory)
# set_target_properties(poisson_DMStag_2D PROPERTIES
#     SOURCES src/poisson/poisson_DMStag_2D.cpp
# )

# Heat equation executable
add_executable(heat_DMStag_2D src/heat/heat_DMStag_2D.cpp)
target_include_directories(heat_DMStag_2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(heat_DMStag_2D PRIVATE PkgConfig::PETSC MPI::MPI_CXX)

# Add tests subdirectory
add_subdirectory(tests/unit/poisson)
add_subdirectory(tests/unit/heat)

